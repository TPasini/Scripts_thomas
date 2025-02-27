import h5py
import argparse

def close_h5_file(filename):
    try:
        with h5py.File(filename, 'a') as h5_file:
            h5_file.close()
            print(f"File {filename} has been closed.")
    except Exception as e:
        print(f"An error occurred while closing {filename}: {e}")

def main():
    parser = argparse.ArgumentParser(description='Force closing of h5 file.')
    parser.add_argument('-h5', '--h5', dest='h5', action='store', default='', type=str, help='H5parm to close.')

    args = parser.parse_args()
    h5file = args.h5

    close_h5_file(h5file)

if __name__ == "__main__":
    main()
