import re
import os
from pathlib import Path

# Calculate runtime
def read_runtimes(result_dir):
    runtimes = []

    pattern = r'took (\d+) seconds'
    with open(result_dir, "r") as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                runtimes.append(int(match.group(1)))
    return sum(runtimes)

# Disk usage used
def get_disk_usage(path):
    total = 0
    for rowsum_dir in Path(path).iterdir():
        if rowsum_dir.is_dir():
            size = sum(f.stat().st_size for f in rowsum_dir.rglob("*") if f.is_file())
            size_MB = size / (1024**2)
            total += size_MB
    return total

# Total QTS count without equivalences
def total_QTS_count(result_dir):
    pattern = r'The function found a total of (\d+) sequences'
    with open(result_dir, 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                return int(match.group(1))
    print('ERROR: Total QTS before equivalence not found')
    exit()
    
# Total QTS count after equivalences
def reduced_QTS_count(result_dir):
    pattern = r'count after equivalences (\d+)'
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
def create_table(start, end, total, S_equ, M_equ, time, pairs, disk, latex):
    if latex is False:
        width = 13
        with open('results.tab', 'w') as file:
            file.write('n'.ljust(width, ' '))
            file.write('Total (s)'.ljust(width, ' '))
            file.write('S_{equ}'.ljust(width, ' '))
            file.write('M_{equ}'.ljust(width, ' '))
            file.write('Time'.ljust(width, ' '))
            file.write('Pairs'.ljust(width, ' '))
            file.write('Disk usage (MB)'.ljust(width, ' '))
            file.write('\n')
            
            for i, n in enumerate(range(start, end)):
                file.write(str(i+1).ljust(width, ' '))
                file.write(str(total[i]).ljust(width, ' '))
                file.write(str(S_equ[i]).ljust(width, ' '))
                file.write(str(M_equ[i]).ljust(width, ' '))
                file.write(str(time[i]).ljust(width, ' '))
                file.write(str(pairs[i]).ljust(width, ' '))
                if disk[i] < 10:
                    file.write(str(round(disk[i], 1)).ljust(width, ' ') + '\n')
                else:
                    file.write(str(round(disk[i])).ljust(width, ' ') + '\n')
    else:
        with open('results.tab', 'w') as file:
            file.write(f'$n$ & Total & $S_{{\\text{{equ}}}}$ & $M_{{\\text{{equ}}}}$ & Time (s) & Pairs & Disk space (MB)\\\\\n')
            for i, n in enumerate(range(start, end)):
                if disk[i] < 10:
                    file.write(f'{i+1} & {total[i]} & {S_equ[i]} & {M_equ[i]} & {time[i]} & {pairs[i]} & {round(disk[i], 1)}\\\\\n')
                else: 
                    file.write(f'{i+1} & {total[i]} & {S_equ[i]} & {M_equ[i]} & {time[i]} & {pairs[i]} & {round(disk[i])}\\\\\n')



# Start execution

start = input("This script collects data from a computation generate given lengths of QTS.\nStart: ")
end = input("End: ")
seqtype = input("Sequence type (qts/wts/ws): ")
latex = True if input("Tab-separated values (t) or LaTeX formatting (l)? ") == 'l' else False

runtime=[]
disk_usage=[]
QTS_total=[]
QTS_reduced=[]
QTS_hadamard_reduced=[]
pairs=[]

for i, n in enumerate(range(int(start), int(end)+1)):
    filePath = "./results/pairs/" + seqtype + "/find_" + str(n)
    result_dir = filePath + "/result.log"

    runtime.append(read_runtimes(result_dir))
    disk_usage.append(get_disk_usage(filePath))
    QTS_total.append(total_QTS_count(result_dir))
    QTS_reduced.append(reduced_QTS_count(result_dir))
    QTS_hadamard_reduced.append(hadamard_reduced_QTS_count(filePath))
    pairs.append(count_pairs(seqtype, n))
    
    print(f'=========================== Length {n} ===========================')
    print(f'Runtime: {runtime[i]} seconds')
    print(f'Disk usage: {disk_usage[i]} MB')
    print(f'QTS without equivalence: {QTS_total[i]}')
    print(f'QTS after sequence equivalence: {QTS_reduced[i]}')
    print(f'QTS after Hadamard equivalence: {QTS_hadamard_reduced[i]}')
    print(f'Total pairs generated: {pairs[i]}\n')

create_table(int(start), int(end), QTS_total, QTS_reduced, QTS_hadamard_reduced, runtime, pairs, disk_usage, latex)
    
