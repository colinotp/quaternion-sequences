import re

n = input("This script collects the total runtime taken to generate a given length of QTS/PQS.\nLength: ")

filePath = "./results/pairs/wts/find_"
directory = filePath + n + "/result.log"


runtimes = []
pattern = r'took (\d+) seconds'
with open(directory, "r") as file:
    for line in file:
        match = re.search(pattern, line)
        if match:
            runtimes.append(int(match.group(1)))

print("Runtimes:")
print(runtimes)
total = sum(runtimes)
print(f"Total runtime: {total} seconds")