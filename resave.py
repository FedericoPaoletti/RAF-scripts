f=open("network_output1.txt","r")
lines=f.readlines()
f.close()

f=open("network_output1.txt","w")
for line in lines:
	f.write(line)
	if line.startswith("//STOP HERE//"):
		break
f.close()
