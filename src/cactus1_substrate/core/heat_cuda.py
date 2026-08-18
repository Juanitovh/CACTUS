"""Optional CUDA implementation of the heat propagation in `workers.meta_grid`.

Import is deliberately lazy and guarded: nothing here is touched unless a caller asks for the GPU path, so
a CPU-only or CUDA-less install behaves exactly as before.

## Why this can be parallelised at all

`meta_grid.propagate_heat_iterations` looks sequential -- it sweeps slices and writes slice `i-1` inside
the loop over `i`. It is not. The write to slice `i-1` happens *after* the read of slices `i-1, i, i+1`, and
no other step writes slice `i-1`, so every read in an iteration sees the pre-iteration array. That makes it
a Jacobi update, and each voxel is independent.

## Quirks that must be reproduced exactly

Verified bit-for-bit against the CPU sweep; each of these was necessary to get there.

1. Per iteration only slices `1 .. lenX-3` are written. The sweep runs `i` to `lenX-2` but writes `i-1`,
   so the last computed slice is discarded.
2. Slice 0 is written with zeros -- it receives the initial `last_slice`, which is never assigned a
   computed value.
3. Every written slice is zero on its own `j` and `k` borders: `propagate_all_heat_2d` builds a fresh
   zero array and its loops skip the borders.
4. Slices `lenX-2` and `lenX-1` are never written by the sweep, so they carry their previous values
   *through the iterations* -- and are read by the stencil at slice `lenX-3`.
5. After all iterations, `grid[0:kernel_size//2] = 0` and `grid[-kernel_size//2:] = 0`. Note that
   `-kernel_size//2` is `(-3)//2 == -2`, not `-1`: unary minus binds looser than `//`, and floor division
   rounds toward negative infinity. So the front loses ONE slice and the back loses TWO. That asymmetry is
   almost certainly unintended upstream, but it is the behaviour, and it is why slice `lenX-2` ends up zero
   despite never being written.

## The stencil

`handle_temperature` builds a dict of neighbour-type counts. The dict is used for exactly two things: the
largest nonzero neighbour type id (`my_max` over its keys -- the largest id, NOT the modal type) and the
count of that one type. Both come out of fixed passes over the 27 neighbours, so the device version needs
no dict and no allocation. Equivalence was property-tested against the dict version on random grids
including sparse, single-type and tie-heavy cases: 0 mismatches.
"""
from __future__ import annotations

import numpy as np


def is_available():
    """Can we actually run on a GPU here? Never raises."""
    try:
        from numba import cuda
        return bool(cuda.is_available())
    except Exception:
        return False


def _build():
    """Compile the kernels on first use, so importing this module stays cheap."""
    from numba import cuda

    @cuda.jit(device=True, inline=True)
    def _stencil(T, C, i, j, k):
        current_type = C[i, j, k]

        predominant = 0.0
        any_nonzero = False
        for idx in range(27):
            it = C[i + idx % 3 - 1, j + idx // 3 % 3 - 1, k + idx // 9 % 3 - 1]
            if it != 0:
                any_nonzero = True
                if it > predominant:
                    predominant = it
        if not any_nonzero:
            return 0.0, 0.0

        count_pre = 0
        for idx in range(27):
            if C[i + idx % 3 - 1, j + idx // 3 % 3 - 1, k + idx // 9 % 3 - 1] == predominant:
                count_pre += 1

        if current_type == 0:
            dominant = predominant
        else:
            dominant = predominant if count_pre > 3 else current_type

        suma = 0.0
        for idx in range(27):
            di, dj, dk = idx % 3 - 1, idx // 3 % 3 - 1, idx // 9 % 3 - 1
            if C[i + di, j + dj, k + dk] == dominant:
                suma += T[i + di, j + dj, k + dk]
        return suma / 27.0, dominant

    @cuda.jit
    def _heat_kernel(T_in, C_in, T_out, C_out, hi_slice):
        i, j, k = cuda.grid(3)
        lenY, lenZ = T_in.shape[1], T_in.shape[2]
        if i < 1 or i > hi_slice:
            return
        if j < 1 or j >= lenY - 1 or k < 1 or k >= lenZ - 1:
            return
        t, c = _stencil(T_in, C_in, i, j, k)
        T_out[i, j, k] = t
        C_out[i, j, k] = c

    return _heat_kernel


_KERNEL = None


def fits_on_device(grid_temperature, grid_class, margin=1.15):
    """Will the four device arrays fit in free GPU memory? Never raises.

    The kernel double-buffers, so it needs 2x the two grids, plus headroom. A caller that skips this
    check gets a CUDA out-of-memory error instead of a clean fall back to the CPU sweep, which on a
    smaller card is the common case rather than an edge one.
    """
    try:
        from numba import cuda
        need = int((grid_temperature.nbytes + grid_class.nbytes) * 2 * margin)
        free, _total = cuda.current_context().get_memory_info()
        return need <= free, need, free
    except Exception:
        return False, 0, 0


def propagate_heat_iterations_cuda(maxIter, grid_temperature, grid_class, kernel_size=3):
    """Drop-in replacement for `meta_grid.propagate_heat_iterations`, on the GPU.

    Same signature minus the stencil function argument (the device stencil is fixed), same in-place
    contract, and bit-identical results. Raises if CUDA is unavailable -- callers should check
    `is_available()` and fall back.
    """
    global _KERNEL
    from numba import cuda

    if _KERNEL is None:
        _KERNEL = _build()

    lenX = grid_temperature.shape[0]
    hi_slice = lenX - 2 - kernel_size // 2      # last slice the CPU sweep actually writes

    # CUDA caps blocks per grid at 65535 in y and z. With the block shape below that allows 2**19 voxels
    # per axis, far beyond anything this code produces, but fail clearly rather than silently truncating.
    tpb = (4, 8, 8)
    bpg = tuple((n + t - 1) // t for n, t in zip(grid_temperature.shape, tpb))
    if bpg[1] > 65535 or bpg[2] > 65535:
        raise ValueError(
            f"grid {grid_temperature.shape} needs {bpg} blocks; CUDA allows at most 65535 in y and z")

    d_T_in = cuda.to_device(np.ascontiguousarray(grid_temperature))
    d_C_in = cuda.to_device(np.ascontiguousarray(grid_class))
    d_T_out = cuda.device_array_like(d_T_in)
    d_C_out = cuda.device_array_like(d_C_in)

    for _ in range(maxIter):
        # zero-filled output reproduces quirks 2 and 3 without a separate pass
        d_T_out[:] = 0
        d_C_out[:] = 0
        _KERNEL[bpg, tpb](d_T_in, d_C_in, d_T_out, d_C_out, hi_slice)
        # quirk 4: the two trailing slices are never written, so carry them forward
        d_T_out[lenX - 2:] = d_T_in[lenX - 2:]
        d_C_out[lenX - 2:] = d_C_in[lenX - 2:]
        d_T_in, d_T_out = d_T_out, d_T_in
        d_C_in, d_C_out = d_C_out, d_C_in
    cuda.synchronize()

    T = d_T_in.copy_to_host()
    C = d_C_in.copy_to_host()
    # quirk 5, using the upstream expressions verbatim
    T[0:kernel_size // 2] = 0
    T[-kernel_size // 2:] = 0
    C[0:kernel_size // 2] = 0
    C[-kernel_size // 2:] = 0
    grid_temperature[...] = T
    grid_class[...] = C
    return grid_temperature, grid_class
