import re
from pathlib import Path

def read_runtimes(directory):
    runtimes = []

    pattern = r'took (\d+) seconds'
    with open(directory, "r") as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                runtimes.append(int(match.group(1)))
    return runtimes

def get_disk_usage(path):
    total = 0
    for rowsum_dir in Path(path).iterdir():
        if rowsum_dir.is_dir():
            size = sum(f.stat().st_size for f in rowsum_dir.rglob("*") if f.is_file())
            size_MB = size / (1024**2)
            total += size_MB
            print(f"Size of {rowsum_dir}: {size_MB}MB\n")
    print(f"Total disk usage for size {n}: {total}MB\n")
    return total

def create_table():
    print("TODO: Finish this")


start = input("This script collects data from a computation generate given lengths of QTS/PQS.\nStart: ")
end = input("End: ")

for n in range(int(start), int(end)+1):
    filePath = "./results/pairs/wts/find_" + str(n)
    result_dir = filePath + "/result.log"

    runtimes = read_runtimes(result_dir)
    disk_usage = get_disk_usage(filePath)
