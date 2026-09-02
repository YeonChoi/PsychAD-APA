import sys
from collections import Counter


def main(id):
    d_cnter = Counter()
    for line in sys.stdin:
        l_fields = line.strip().split()

        for field in reversed(l_fields):
            if field.startswith('XS:Z'):
                d_cnter[field] += 1
                break

    print(id, d_cnter, sep='\t')

if __name__=='__main__':
    id = sys.argv[1]
    main(id)

