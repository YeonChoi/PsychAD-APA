import subprocess
import sys
import pandas as pd

def read_barcode_list(bc_path):
    df = pd.read_table(bc_path, index_col=False, header=None)
    if len(df)!=len(set(df[0])):
        raise ValueError(f'Barcodes should be unique! counts are: total: {len(df)}, unique: {set(df[0])}')
    return dict(zip(df[0], df[1]))

def get_CB(items):
    for item in reversed(items[11:]):
        if item.startswith('CB:Z'):
            return item
    return None

def write_splitted_sam(file_objs, cb_to_celltypes):
    for line in sys.stdin:
        if line.startswith('@'):
            sys.stdout.write(line)
            continue
        items = line.strip().split('\t')
        cb = get_CB(items)
        if cb is None: continue

        try:
            cell_type = cb_to_celltypes[cb]
        except KeyError:
            raise KeyError(f'cell barcode {cb} does not exists in cell type file!')

        print(line, file=file_objs[cell_type], end='')

    for fi in file_objs.values():
        fi.close()

def open_output_sam_files(cb_to_celltypes, out_prefix):
    cell_types = set(cb_to_celltypes.values())
    opened_file_objs = {}

    for ct in cell_types:
        out_path = out_prefix + ct + '.sam'
        ct_sam_file = open(out_path, 'a')
        opened_file_objs[ct] = ct_sam_file

    return opened_file_objs


if __name__ == "__main__":
    bc_path = sys.argv[1]
    out_prefix = sys.argv[2]

    cb_to_celltypes = read_barcode_list(bc_path)
    opened_file_objs = open_output_sam_files(cb_to_celltypes, out_prefix)
    write_splitted_sam(opened_file_objs, cb_to_celltypes)
