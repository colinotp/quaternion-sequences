import time, sys, os.path
from subprocess32 import call

decomps_len = [0, 1, 1, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 3, 1, 2, 2, 2, 3, 3, 2, 3, 2, 2, 1, 4, 3, 4, 3, 2, 2, 4, 1, 4, 4, 3, 4, 5, 3, 4, 2, 3, 4, 5, 2, 6, 2, 4, 2, 6, 5, 4, 5, 4, 5, 6, 1, 6, 5, 4, 3, 6, 3, 8, 2, 5, 6, 7, 4, 6, 5]

verbose = True
uselogs = False

if not uselogs:
	logfile = open("/dev/null", "w")

if len(sys.argv) < 2:
	print "Need order of search to run, and optionally a case # to solve and the compression factor"
	print "Passing -1 as the case # will solve all cases"
	quit()

n = int(sys.argv[1])

casetosolve = int(sys.argv[2]) if len(sys.argv) > 2 else -1

def runstr(n, c, k):
	if k == -1 and c == -1:
		return "%s" % n
	elif k == -1:
		return "%s.%s" % (n, c)
	else:
		return "%s.%s.%s" % (n, c, k)

inname = "input/compstring/%s.in"
logname = "output/compstring/%s.log"
exhaustname = "exhaust/%s"
timingsname = "timings/%s.solvetime"

if not os.path.exists("exhaust"): call(["mkdir", "exhaust"])
if not os.path.exists("output"): call(["mkdir", "output"])
if not os.path.exists("output/compstring"): call(["mkdir", "output/compstring"])
if not os.path.exists("results"): call(["mkdir", "results"])
if not os.path.exists("results/compstring"): call(["mkdir", "results/compstring"])
if not os.path.exists("timings"): call(["mkdir", "timings"])

if len(sys.argv) > 3:
	d = int(sys.argv[3])
else:
	d = 1

l = n/d

print "ORDER {0}".format(n)

cases = range(decomps_len[n])

for c in cases:

	if casetosolve != -1 and casetosolve != c:
		continue

	if os.path.exists(exhaustname % runstr(n, c, -1)):
		call("mv {0} {1}".format(exhaustname % runstr(n, c, -1), exhaustname % runstr(n, c, -1) + ".orig").split(" "))
	call("touch {0}".format(exhaustname % runstr(n, c, -1)).split(" "))

	if not os.path.exists("matchedseqns/{0}.{1}.{2}.inequiv".format(n, c, l)):
		print "Need file matchedseqns/{0}.{1}.{2}.inequiv".format(n, c, l)
		continue

	f = open("matchedseqns/{0}.{1}.{2}.inequiv".format(n, c, l), "r")
	lines = f.readlines()
	f.close()

	start = time.time()
	totalnumsols = 0
	totalfilttime = 0
	totalsolved = 0
	totalsat = 0
	totaltime = 0
	maxtime = 0

	for li in range(len(lines)):

		pos = lines[li][:-2].split(" ")
		pos = map(int, pos)
		assert(len(pos)==4*l)

		compstring = "{0},".format(d)
		for x in pos:
			compstring += "{0},".format(x)
		compstring = compstring[:-1]

		command = "./maplesat_static -no-pre -order={0} ".format(n) + inname % runstr(n, -1, -1) + " -exhaustive=" + exhaustname % runstr(n, c, -1) #+ " " + resultname % runstr(n, c, k)
		command += " -filtering"
		command += " -compstring=" + compstring

		if verbose:
			print command

		if uselogs:
			thisstart = time.time()
			logfile = open(logname % runstr(n, c, li), "w")
			retval = call(command.split(" "), stdout=logfile)
			logfile.close()
			thistime = time.time()-thisstart
			totaltime += thistime
			if thistime > maxtime:
				maxtime = thistime

			assert(retval == 20 or retval == 21)

			numsols = 0
			filttime = 0
			logfile = open(logname % runstr(n, c, li), "r")
			for line in logfile.readlines():
				if "NUMSOLS" in line:
					numsols = int(line.split(" ")[-1])
					if verbose and numsols > 0:
						print "NUMSOLS: {0}".format(numsols)
				if "filtering   checks:" in line:
					filttime = float(line.split(" ")[-3])
					if verbose:
						print "filttime: {0}".format(filttime)
			logfile.close()
			totalfilttime += filttime
			if numsols > 0:
				totalnumsols += numsols
				totalsat += 1
		else:
			call(command.split(" "), stdout=logfile)

		totalsolved += 1

	if uselogs:
		print "  Case {0}: {1}/{2} SAT instances with {3} total solutions solved in %.2f seconds (max: %.2f, avg: %.2f) realtime: %.2f seconds filttime: %.2f seconds".format(c, totalsat, totalsolved, totalnumsols) % (totaltime, maxtime, totaltime/totalsolved if totalsolved != 0 else 0, time.time()-start, totalfilttime)
	else:
		print "  Case {0}: {1} SAT instances solved in %.2f seconds".format(c, totalsolved) % (time.time()-start)

	f = open(timingsname % runstr(n, c, -1), "w")
	if uselogs:
		f.write("%.2f\n" % totaltime)
	else:
		f.write("%.2f\n" % (time.time()-start))
	f.close()
