import argparse
import sys
import os

def parser():
    parser = argparse.ArgumentParser(description="Generates some sample documents in a dataset")
    parser.add_argument('--output_path', required=True, action="store", type=str, help="Path to output file")
    return parser.parse_args()


if __name__=='__main__':

        args = parser()

        outdir = args.output_path

        if not os.path.exists(args.output_path):
            os.makedirs(args.output_path)
        for name in ['one', 'two', 'three']:
            with open(os.path.join(outdir, name + ".txt"), 'w') as out:
                out.write("%s, 1.3\n" % name)
