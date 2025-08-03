#!/usr/bin/env python3
#
# Copyright (c) 2025 Mate Soos
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import argparse
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--verbose", help="Verbose debug printing", action="store_const",
        const=True)
    parser.add_argument("inputFile", help="input File (in Weighted CNF format)")
    parser.add_argument("outputFile", help="output File (in Weighted CNF format)")
    args = parser.parse_args()

    with open(args.outputFile, 'w') as out:
      with open(args.inputFile, 'r') as f:
        for line in f:
          line = line.strip()
          if line.startswith("c p weight"):
            line = re.sub(r'\s+', ' ', line)
            if not re.fullmatch(r'c p weight \S+ \S+ \S+ 0', line):
                print(f"Error: Invalid line format: {line}")
                print(f"Expected format: 'c p weight <var1> <real weight> <complex weight> 0'")
                print(f"Instead got: {line}")
                exit(1)
            line = re.sub(r'c p weight (\S+) (\S+) (\S+) 0', r'c p weight \1 \2 + \3i 0', line)
            out.write(line + "\n")
          else:
            out.write(line + "\n")
