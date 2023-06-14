import sys, os.path
from subprocess32 import call

if len(sys.argv) < 2:
	print "need order (or min and max order) of cases to generate"
	quit()

n1 = int(sys.argv[1])
if len(sys.argv) > 2:
	n2 = int(sys.argv[2])
else:
	n2 = n1

inname = "input/compstring/%s.in"

if not os.path.exists("input"): call(["mkdir", "input"])
if not os.path.exists("input/compstring"): call(["mkdir", "input/compstring"])

def runstr(n, c, k):
	if k == -1 and c == -1:
		return "%s" % n
	elif k == -1:
		return "%s.%s" % (n, c)
	else:
		return "%s.%s.%s" % (n, c, k)

for n in range(n1, n2+1):
	#indices = [0*n+1, 1*n+1, 2*n+1, 3*n+1]
	f = open(inname % runstr(n, -1, -1), "w")
	f.write("p cnf {0} 1\n".format(4*n))
	f.write("{0} -{0} 0\n".format(4*n))

	f.close()
