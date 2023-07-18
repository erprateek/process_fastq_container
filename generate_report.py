import os
import sys
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",'--fastqc_outputs_dir')
    parser.add_argument("-t"."--trimgalore_output_dir")
    parser.add_argument("-o","--output_file")
    
    

if __name__ == '__main__':
    main()