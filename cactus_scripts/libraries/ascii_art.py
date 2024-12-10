


                                                                               

cactus_ascii = """
 ░▒▓██████▓▒░ ░▒▓██████▓▒░ ░▒▓██████▓▒░▒▓████████▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓███████▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        
░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        
░▒▓█▓▒░      ░▒▓████████▓▒░▒▓█▓▒░        ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░░▒▓██████▓▒░  
░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░ 
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░ 
 ░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░░▒▓██████▓▒░  ░▒▓█▓▒░    ░▒▓██████▓▒░░▒▓███████▓▒░  
"""


def display_title():
    """Print ASCII art for the program."""
    print(cactus_ascii)
                                                                               
                                                                               

                                                                               
def print_message(message):
    # Determine the longest line length
    max_line_length = max(len(line) for line in message.split('\n'))  
    padding = 2  # Padding around the text
    box_width = max_line_length + padding * 2

    # Bold box characters
    top_left = "┏"
    top_right = "┓"
    bottom_left = "┗"
    bottom_right = "┛"
    horizontal = "━"
    vertical = "┃"

    # Top border
    print(f"{top_left}{horizontal * box_width}{top_right}")
    
    # Message lines with padding
    for line in message.split('\n'):
        padded_line = line.center(max_line_length)  # Center the text
        print(f"{vertical} {' ' * padding}{padded_line}{' ' * (padding-1)}{vertical}")
    
    # Bottom border
    print(f"{bottom_left}{horizontal * box_width}{bottom_right}")

# Example Usage
                                                                               

def print_important_message(message):
    # Determine the longest line length
    max_line_length = max(len(line) for line in message.split('\n'))
    #round it to odd number
    if max_line_length % 2 == 1:
        max_line_length += 1
    padding = 4  # Extra padding for important messages
    box_width = max_line_length + padding * 2

    # Bold box characters
    top_left = "╔"
    top_right = "╗"
    bottom_left = "╚"
    bottom_right = "╝"
    horizontal = "═"
    vertical = "║"

    # Top border with "IMPORTANT" tag
    print(f"{top_left}{horizontal * ((box_width - 9) // 2)} IMPORTANT {horizontal * ((box_width - 9) // 2)}{top_right}")

    # Message lines with padding
    for line in message.split('\n'):
        padded_line = line.center(max_line_length)  # Center the text
        print(f"{vertical} {' ' * padding}{padded_line}{' ' * padding}{vertical}")

    # Bottom border
    print(f"{bottom_left}{horizontal * (box_width +1)}{bottom_right}")


def print_error_message(message):
    # Determine the longest line length
    max_line_length = max(len(line) for line in message.split('\n'))
    if max_line_length % 2 == 1:
        max_line_length += 1
    padding = 6  # Extra padding for dramatic effect
    box_width = max_line_length + padding * 2

    # Box characters
    top_left = "╔"
    top_right = "╗"
    bottom_left = "╚"
    bottom_right = "╝"
    horizontal = "═"
    vertical = "║"

    # Over-the-top top banner
    print("\n" + "✖" * (box_width + 6))  # Top dramatic border
    print(f"!{' ERROR '.center(box_width + 2, '!')}!")

    # Top border
    print(f"{top_left}{horizontal * (box_width+2)}{top_right}")

    # Message lines with padding
    for line in message.split('\n'):
        padded_line = line.center(max_line_length)  # Center the text
        print(f"{vertical} {' ' * padding}{padded_line}{' ' * padding} {vertical}")

    # Bottom border
    print(f"{bottom_left}{horizontal * (box_width+2)}{bottom_right}")

    # Over-the-top bottom banner
    print("✖" * (box_width + 6))

# Example Usage
# Example Usage

if __name__ == "__main__":
    print_message("This is a simple message \nwith multiple lines.")
    print_important_message("CRITICAL ERROR!\nSystem will shut down in 30 seconds.\nPlease save your work immediately.")
    print_error_message("CRITICA ERROR!\nSystem crash imminent.\nBackup your data and restart immediatelyu!")
                                                                               
                                                                               
                                                                               
