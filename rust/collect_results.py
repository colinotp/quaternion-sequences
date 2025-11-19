#!/usr/bin/env python3

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
def reduced_QTS_count(result_dir, seqtype):
    pattern = r'Found (\d+) ' + seqtype + ' after reducing to equivalence'
    with open(result_dir, 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                return int(match.group(1))
    return -1
    
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
    pattern = r'Generated (\d+) total pairs'
    result_log = "./results/pairs/" + seqtype + "/find_" + str(n) + "/result.log"
    if os.path.isfile(result_log):
        with open(result_log, 'r') as file:
            for line in file:
                match = re.search(pattern, line)
                if match:
                    return int(match.group(1))
    return int(os.popen('./countpairs.sh ' + seqtype + ' ' + str(n)).read())

# Create table from data. Each arg other than n should be a list of length n
def create_table(start, end, Q_equ, W_equ, M_equ, time, qt_equiv_time, hm_equiv_time, pairs, disk, latex, seqtype):
    if latex is False:
        width = 16
    
        print('n'.rjust(2, ' '), end='')
        if seqtype == "wts":
            print('W_{equ}'.rjust(width, ' '), end='')
        print('Q_{equ}'.rjust(width, ' '), end='')
        print('M_{equ}'.rjust(width, ' '), end='')
        print('Time (s)'.rjust(width, ' '), end='')
        print('QT Equ Time (s)'.rjust(width, ' '), end='')
        print('HM Equ Time (s)'.rjust(width, ' '), end='')
        print('Pairs'.rjust(width, ' '), end='')
        print('Disk usage (MB)'.rjust(width, ' '))
        
        for i, n in enumerate(range(start, end+1)):
            print(str(n).rjust(2, ' '), end='')
            if seqtype == "wts":
                print(str(W_equ[i]).rjust(width, ' '), end='')
            print(str(Q_equ[i]).rjust(width, ' '), end='')
            print(str(M_equ[i]).rjust(width, ' '), end='')
            print(str(f"{time[i]:.2f}").rjust(width, ' '), end='')
            print(str(f"{qt_equiv_time[i]:.2f}").rjust(width, ' '), end='')
            print(str(f"{hm_equiv_time[i]:.2f}").rjust(width, ' '), end='')
            print(str(pairs[i]).rjust(width, ' '), end='')
            print(f"{disk[i]:.1f}".rjust(width, ' '))
    else:
        if seqtype == "wts":
            print(f'$n$ & $W_{{\\text{{equ}}}}$ & $Q_{{\\text{{equ}}}}$ & $M_{{\\text{{equ}}}}$ & Time (s) & Pairs & Space \\\\')
        else:
            print(f'$n$ & $Q_{{\\text{{equ}}}}$ & $M_{{\\text{{equ}}}}$ & Time (s) & Pairs & Space \\\\')
        print('\\hline')
        for i, n in enumerate(range(start, end+1)):
            if seqtype == "wts":
                print(f'{n} & {W_equ[i]} & {Q_equ[i]} & {M_equ[i]} & {time[i]:.2f} & {pairs[i]} & {disk[i]:.1f} \\\\')
            else:
                print(f'{n} & {Q_equ[i]} & {M_equ[i]} & {time[i]:.2f} & {pairs[i]} & {disk[i]:.1f} \\\\')



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
WTS_reduced=[]
QTS_hadamard_reduced=[]
pairs=[]

for i, n in enumerate(range(int(start), int(end)+1)):
    filePath = "./results/pairs/" + seqtype + "/find_" + str(n)
    result_dir = filePath + "/result.log"

    if not os.path.isfile(result_dir):
        if verbose:
            print(f"File {result_dir} doesn't exist")
        runtime.append(-1)
        qt_equivalence_time.append(-1)
        hm_equivalence_time.append(-1)
        disk_usage.append(-1)
        QTS_reduced.append(-1)
        WTS_reduced.append(-1)
        QTS_hadamard_reduced.append(-1)
        pairs.append(-1)
        continue

    runtime.append(read_runtimes(result_dir))
    qt_equivalence_time.append(get_qt_equivalence_time(result_dir))
    hm_equivalence_time.append(get_hm_equivalence_time(result_dir))
    disk_usage.append(get_disk_usage(filePath))
    QTS_reduced.append(reduced_QTS_count(result_dir, "qts"))
    WTS_reduced.append(reduced_QTS_count(result_dir, "wts"))
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

if verbose:
    if seqtype == "wts":
        print(f"A total of {sum(WTS_reduced)} sets of QT sequences found up to Williamson-type equivalence")
    print(f"A total of {sum(QTS_reduced)} sets of QT sequences found up to QT equivalence")
    print(f"A total of {sum(QTS_hadamard_reduced)} sets of QT sequences found up to Hadamard equivalence")
    print("")

create_table(int(start), int(end), QTS_reduced, WTS_reduced, QTS_hadamard_reduced, runtime, qt_equivalence_time, hm_equivalence_time, pairs, disk_usage, latex, seqtype)
