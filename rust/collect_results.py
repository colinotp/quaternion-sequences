import re
import os
import sys
from pathlib import Path

# Calculate runtime
def read_runtimes(result_dir):
    runtimes = []

    pattern = r': (\d+\.\d\d) seconds'
    with open(result_dir, "r") as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                runtimes.append(float(match.group(1)))
    return sum(runtimes)
    

# Get time taken to reduce to QT equivalence
def get_qt_equivalence_time(result_dir):
    pattern = r'Reducing to equivalence took (\d+(?:\.\d+)?) seconds'
    with open(result_dir, "r") as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                return float(match.group(1))
        return -1
    
# Get time taken to reduce to Hadamard equivalence
def get_hm_equivalence_time(result_dir):
    pattern = r'Converting to matrices up to Hadamard equivalence took (\d+(?:\.\d+)?) seconds'
    with open(result_dir, "r") as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                return float(match.group(1))
        return -1

# Disk usage used
def get_disk_usage(path):
    total = 0
    for rowsum_dir in Path(path).iterdir():
        if rowsum_dir.is_dir():
            size = sum(f.stat().st_size for f in rowsum_dir.rglob("*") if f.is_file())
            size_MB = size / (1024**2)
            total += size_MB
    return total

# Total QTS count after equivalences
def reduced_QTS_count(result_dir):
    pattern = r'Found (\d+) (\w+) after reducing to equivalence'
    with open(result_dir, 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                return int(match.group(1))
    print('ERROR: QTS after equivalence not found')
    exit()
    
# Total Hadamard matrix count up to hadamard equivalence
def hadamard_reduced_QTS_count(path):
    result_mat = path + '/result.mat'
    try:
        with open(result_mat, 'r') as file:
            count = len(file.readlines())
    except:
        return -1
    return count
    
# Count generated pairs
def count_pairs(seqtype, n):
    return int(os.popen('./countpairs.sh ' + seqtype + ' ' + str(n)).read())

# Create table from data. Each arg other than n should be a list of length n
def create_table(start, end, S_equ, M_equ, time, qt_equiv_time, hm_equiv_time, pairs, disk, latex):
    if latex is False:
        width = 16
    
        print('n'.ljust(width, ' '), end='')
        print('S_{equ}'.ljust(width, ' '), end='')
        print('M_{equ}'.ljust(width, ' '), end='')
        print('Time (s)'.ljust(width, ' '), end='')
        print('QT Equ Time (s)'.ljust(width, ' '), end='')
        print('HM Equ Time (s)'.ljust(width, ' '), end='')
        print('Pairs'.ljust(width, ' '), end='')
        print('Disk usage (MB)'.ljust(width, ' '))
        
        for i, n in enumerate(range(start, end+1)):
            print(str(n).ljust(width, ' '), end='')
            print(str(S_equ[i]).ljust(width, ' '), end='')
            print(str(M_equ[i]).ljust(width, ' '), end='')
            print(str(round(time[i], 2)).ljust(width, ' '), end='')
            print(str(round(qt_equiv_time[i], 2)).ljust(width, ' '), end='')
            print(str(round(hm_equiv_time[i], 2)).ljust(width, ' '), end='')
            print(str(pairs[i]).ljust(width, ' '), end='')
            if disk[i] < 10:
                print(str(round(disk[i], 1)).ljust(width, ' '))
            else:
                print(str(round(disk[i])).ljust(width, ' '))
    else:
        with open('results.tab', 'w') as file:
            print(f'$n$ & $S_{{\\text{{equ}}}}$ & $M_{{\\text{{equ}}}}$ & Time (s) & QT Equ Time (s) & HM Equ Time (s) & Pairs & Disk space (MB)\\\\')
            for i, n in enumerate(range(start, end+1)):
                if disk[i] < 10:
                    print(f'{n} & {S_equ[i]} & {M_equ[i]} & {round(time[i], 2)} & {round(qt_equiv_time[i], 2)} & {round(hm_equiv_time[i], 2)} & {pairs[i]} & {round(disk[i], 1)}\\\\')
                else: 
                    print(f'{n} & {S_equ[i]} & {M_equ[i]} & {round(time[i], 2)} & {round(qt_equiv_time[i], 2)} & {round(hm_equiv_time[i], 2)} & {pairs[i]} & {round(disk[i])}\\\\')



# Start execution
verbose = False
if len(sys.argv) == 6 and sys.argv[5] == "verbose":
    verbose = True
elif len(sys.argv) != 5:
    print("This script collects data from a computation generate given lengths of QTS.")
    print("To generate table for sequences of lengths a-b (inclusive), use args <a> <b> <sequencetype> <t/l>, where t is used to generate a tsv table, and l is used to generate a latex table.")
    print("If 'verbose' is passed as an additional argument, the output will be verbose.")
    exit(0)

start, end, seqtype = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
latex = True if sys.argv[4] == 'l' else False

runtime=[]
qt_equivalence_time=[]
hm_equivalence_time=[]
disk_usage=[]
QTS_reduced=[]
QTS_hadamard_reduced=[]
pairs=[]

for i, n in enumerate(range(int(start), int(end)+1)):
    filePath = "./results/pairs/" + seqtype + "/find_" + str(n)
    result_dir = filePath + "/result.log"

    runtime.append(read_runtimes(result_dir))
    qt_equivalence_time.append(get_qt_equivalence_time(result_dir))
    hm_equivalence_time.append(get_hm_equivalence_time(result_dir))
    disk_usage.append(get_disk_usage(filePath))
    QTS_reduced.append(reduced_QTS_count(result_dir))
    QTS_hadamard_reduced.append(hadamard_reduced_QTS_count(filePath))
    pairs.append(count_pairs(seqtype, n))
    
    if verbose:
        print(f'=========================== Length {n} ===========================')
        print(f'Runtime: {runtime[i]} seconds')
        print(f'Time to reduce to equivalence: {round(qt_equivalence_time[i])} seconds')
        print(f'Disk usage: {disk_usage[i]} MB')
        print(f'QTS after sequence equivalence: {QTS_reduced[i]}')
        print(f'QTS after Hadamard equivalence: {QTS_hadamard_reduced[i]}')
        print(f'Total pairs generated: {pairs[i]}\n')


create_table(int(start), int(end), QTS_reduced, QTS_hadamard_reduced, runtime, qt_equivalence_time, hm_equivalence_time, pairs, disk_usage, latex)
