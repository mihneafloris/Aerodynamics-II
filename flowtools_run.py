import flowtools

gamma = 1.4

out = flowtools.flowisentropic2(gamma,50.0,'sub')
print(out)

out2 = flowtools.flownormalshock2(gamma,2.5,'mach')
print(out2)

