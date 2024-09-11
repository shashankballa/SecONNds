import re
import sys
import os

def process_conv_block(block):
    #remove line breaks from the block
    block_line = block.replace('\n', '')
    #remove whitespace from the block
    block_line = block_line.replace(' ', '')
    # Regex to match Conv2DWrapper call with tmp* variables
    matches = re.findall(r'Conv2DWrapper\((.*?)tmp(\d+),tmp(\d+),(.*?)\);', block_line, re.DOTALL)

    if matches:
        for match in matches:
            # Extract parts of the Conv2DWrapper call
            part1 = match[0]  # Everything before the first tmp (includes conv_ntt and other args)
            num1 = match[1]   # Number for the first tmp (first tmp argument)
            num2 = match[2]   # Number for the second tmp (second tmp argument)
            part4 = match[3]  # Everything after the second tmp (remaining args)

            # Construct the ConvOnlineHeliks block, ensuring the first tmp argument is included
            conv_heliks = f"""ConvOnlineHeliks(conv_ntt,{part1}tmp{num1},tmp{num2},noise_cts_{num2},noise_pts_{num2},secret_share_vec_{num2},tmp{num2}_pts,{part4});"""

            # # Check if there is already an else block or insert the ConvOnlineCheetah
            # if 'else' not in block:
            #     block = block + f'\n  else {{{conv_cheetah}\n  }}'
            # else:
            #     # Append ConvOnlineCheetah inside the existing else block
            #     block = re.sub(r'else\s*\{', f' {conv_cheetah}', block)

            block_line = "if (!use_seconnds) {\n\t\t" + block_line + "\n\t} else {\n\t\t" + conv_heliks + "\n\t}"

            # print(replace_block)

    return block_line

def process_cpp_code(cpp_code):
    # # Pattern to find if(!use_seconnds) blocks
    # pattern = re.compile(r'if\s*\(!use_seconnds\)\s*\{[^}]*\}', re.DOTALL)

    # Pattern to find Conv2DWrapper calls with any arguments
    pattern = re.compile(r'Conv2DWrapper\(.*?\);', re.DOTALL)

    # print the number of matches
    print(f"Number of Conv2DWrapper calls found: {len(re.findall(pattern, cpp_code))}")

    # Apply the process_conv_block function to each found block
    updated_code = pattern.sub(lambda m: process_conv_block(m.group(0)), cpp_code)

    return updated_code

def main(input_file):
    # Read the C++ code from the input file
    with open(input_file, 'r') as file:
        cpp_code = file.read()

    # Process the C++ code
    updated_code = process_cpp_code(cpp_code)

    # Generate the output file name by appending "_conv_on" to the input filename (without extension)
    base_name, ext = os.path.splitext(input_file)
    output_file = f"{base_name}_conv_on{ext}"

    # Write the updated C++ code to the output file
    with open(output_file, 'w') as file:
        file.write(updated_code)

    print(f"Processed code has been saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python process_cpp.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]

    main(input_file)
