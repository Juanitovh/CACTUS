"""CLI entry point for cactus1-substrates."""

import argparse
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="cactus1-substrates",
        description="CACTUS — Computational Axonal Configurator for Tailored and Ultradense Substrates",
    )
    sub = parser.add_subparsers(dest="command")

    # --- init ---
    p_init = sub.add_parser("init", help="Initialize fibre placements")
    p_init.add_argument("-config_file", type=str, required=True, help="Configuration file")

    # --- optimize ---
    p_opt = sub.add_parser("optimize", help="Run global joint optimization")
    p_opt.add_argument("-config_file", type=str, required=True, help="Configuration file")
    p_opt.add_argument("-file", type=str, default=None, help="*.init file to process")

    # --- grow ---
    p_grow = sub.add_parser("grow", help="Fibre radial growth and meshing")
    p_grow.add_argument("-config_file", type=str, required=True, help="Configuration file")
    p_grow.add_argument("-substep", choices=["growth", "mesh"], required=True, help="Substep to run")
    p_grow.add_argument("-run_case", choices=["test", "missing", "all"], required=True, help="Run case")
    p_grow.add_argument("-file", type=str, default=None, help="*.init file to process")

    # --- quick-mesh ---
    p_qm = sub.add_parser("quick-mesh", help="Quick mesh generation from strand file")
    p_qm.add_argument("-file", type=str, required=True, help="Input strand file")
    p_qm.add_argument("-file_out", type=str, default="test2.ply", help="Output mesh file")
    p_qm.add_argument("-g_ratio", type=float, default=1.0, help="Growth ratio for scaling radii")
    p_qm.add_argument("-strand_id", type=str, default=None, help="Comma-separated strand IDs")
    p_qm.add_argument("-decimate", type=float, default=0.0, help="Decimation percentage")
    p_qm.add_argument("--parallel", action="store_true", help="Use parallel processing")
    p_qm.add_argument("-n_circles", type=int, default=12, help="Number of circles (parallel mode)")
    p_qm.add_argument("-output", type=str, default="output_mesh.ply", help="Output file (parallel mode)")

    # --- monitor ---
    p_mon = sub.add_parser("monitor", help="Monitor optimization progress")
    p_mon.add_argument("-folder", type=str, required=True, help="Folder containing log files")

    # --- convert ---
    p_conv = sub.add_parser("convert", help="Convert NFG format to CACTUS format")
    p_conv.add_argument("-folder", type=str, required=True, help="Folder with NFG strand files")
    p_conv.add_argument("-outfile", type=str, required=True, help="Output file")

    # --- paste-mesh ---
    p_paste = sub.add_parser("paste-mesh", help="Paste individual strand meshes into a single mesh")
    p_paste.add_argument("-file", type=str, required=True, help="Strand file")
    p_paste.add_argument("-n_erode", type=str, default="0", help="Erosion index")
    p_paste.add_argument("-sim_vol", type=str, default="simulations", help="Folder of pasting")
    p_paste.add_argument("-inn_out", type=str, default="outer", help="Inner or outer")
    p_paste.add_argument("-decimate", type=float, default=0, help="Decimation [0,1)")
    p_paste.add_argument("-patch_flip", type=int, default=1, help="0/1")
    p_paste.add_argument("-colorless", action="store_true", help="Colorless vertex")

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "init":
        from cactus1_substrate.pipeline.initialization import main as run
        # Build a namespace matching what the pipeline expects
        run_args = argparse.Namespace(config_file=args.config_file)
        run(run_args)

    elif args.command == "optimize":
        from cactus1_substrate.pipeline.optimization import main as run
        run_args = argparse.Namespace(config_file=args.config_file, file=args.file)
        run(run_args)

    elif args.command == "grow":
        from cactus1_substrate.pipeline.growth_mesh import main as run
        run_args = argparse.Namespace(
            config_file=args.config_file,
            substep=args.substep,
            run_case=args.run_case,
            file=args.file,
        )
        run(run_args)

    elif args.command == "quick-mesh":
        if args.parallel:
            from cactus1_substrate.tools.quick_mesh_parallel import main as run
            run_args = argparse.Namespace(
                file=args.file,
                output=args.output,
                g_ratio=args.g_ratio,
                n_circles=args.n_circles,
                decimate=args.decimate,
                strand_id=args.strand_id,
            )
            run(run_args)
        else:
            from cactus1_substrate.tools.quick_mesh import main as run
            run_args = argparse.Namespace(
                file=args.file,
                file_out=args.file_out,
                g_ratio=args.g_ratio,
                strand_id=args.strand_id,
                decimate=args.decimate,
            )
            run(run_args)

    elif args.command == "monitor":
        from cactus1_substrate.tools.optim_logger import main as run
        run_args = argparse.Namespace(folder=args.folder)
        run(run_args)

    elif args.command == "convert":
        from cactus1_substrate.tools.nfg2cactus import process_my_strands
        process_my_strands(args.folder, args.outfile)

    elif args.command == "paste-mesh":
        from cactus1_substrate.tools.paste_mesh_dataset import main as run
        run(args)
