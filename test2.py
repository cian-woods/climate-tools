import commands,os,glob,sys

Model          = 'CCSM4'
Dir            = '../rcp85/%s/day/plev/va' % (Model)
Files          = glob.glob('%s/*' % (Dir))
Times          = []
for fname in Files:
	print fname
	bashcommand1   = '''ncdump -v time %s | grep -n time | awk -F  ":" '{print $1}' | tail -1''' % (fname)
	status,output1 = commands.getstatusoutput(bashcommand1)
	bashcommand2   = '''ncdump -v time %s | tail -n +%s''' % (fname,output1)
	status,output2 = commands.getstatusoutput(bashcommand2)
	lines = output2.split('\n')

	L = []
	line    = lines[0].split(',')[:-1]
	line[0] = line[0][7:]
	line    = [float(l) for l in line]
	L       = L + line
	for i in range(1,len(lines)-2,1):
		line = lines[i].split(',')
		L = L + [float(line[l]) for l in range(len(line)-1)]
	line = lines[-2].split(',')
	line[-1] = line[-1][0:-2]
	line = [float(l) for l in line]
	L    = L + line
	Times.append(L)
