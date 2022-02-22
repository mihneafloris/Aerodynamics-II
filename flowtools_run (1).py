import flowtools
import numpy as np
import matplotlib.pyplot as plt

gamma = 1.4

M4 = 18
M5 = 10

#X location of the beginning and ending points of the nozzle block
x_nozzle_begin = 44.8
x_nozzle_end = 194.8
#Location of the 1st throat
x_1st_throat = 65.0
#X location of the beginning and ending points of the flexible diffuser
x_diffuser_begin = 410.0
x_diffuser_end = 760.0

x_shock3 = 630.0
x_shock4 = 530.0
x_shock5 = x_diffuser_begin

#Generate data of geometry of the nozzle block (44.8 mm ≤ x ≤ 194.8 mm)
nb_dat = np.genfromtxt("nozzle_block.txt", comments='%')
nb_AAt = nb_dat[:,2]

#Generate data of geometric property of the flexible diffuser (410mm ≤ x ≤ 760mm)
fd_dat = np.genfromtxt("flexible_diffuser.txt", comments='%')
#Calculate the ratio of the area of the flexible diffuser to the 2nd throat
def DiffuserAreaRatio(x, M4, M5):
    vals = fd_dat[:,2:5][np.where(fd_dat[:,0]==M4)][np.where(fd_dat[:,1][np.where(fd_dat[:,0]==M4)]==M5)].flatten()
    h = vals[1]*x+vals[2]
    AAt2 = h/vals[0]
    return AAt2

#Generate data from test files of part 1.1 and 1.2
tf1_dat = np.genfromtxt("testfile_1.txt", comments='%')
tf2_dat = np.genfromtxt("testfile_2.txt", comments='%')
#Find the beginning and ending indeces of the nozzle
idx_nb_begin = np.argwhere(tf1_dat[:,0]>=x_nozzle_begin).flatten()[0]
idx_nb_end = np.argwhere(tf1_dat[:,0]<=x_nozzle_end).flatten()[-1]
#Location and pressure ratio of the nozzle
nb_mm = tf1_dat[idx_nb_begin:idx_nb_end+1,0]
tf1_ppt = tf1_dat[idx_nb_begin:idx_nb_end+1,1]
tf2_ppt = tf2_dat[idx_nb_begin:idx_nb_end+1,1]

#Generate data from test files of part 2.3, 2.4 and 2.5
tf3_dat = np.genfromtxt("testfile_3.txt", comments='%')
tf4_dat = np.genfromtxt("testfile_4.txt", comments='%')
tf5_dat = np.genfromtxt("testfile_5.txt", comments='%')
#Find the beginning and ending indeces of the diffuser
idx_fd_begin = np.argwhere(tf3_dat[:,0]>=x_diffuser_begin).flatten()[0]
idx_fd_end = np.argwhere(tf3_dat[:,0]<=x_diffuser_end).flatten()[-1]
#Calculate the geometry of the flexible diffuser
fd_mm = tf3_dat[idx_fd_begin:idx_fd_end+1,0]
fd_AAt = DiffuserAreaRatio(fd_mm, M4, M5)
#Location and pressure ratio of the whole wind tunnel
wt_mm = tf3_dat[:,0]
tf3_ppt = tf3_dat[:,1]
tf4_ppt = tf4_dat[:,1]
tf5_ppt = tf5_dat[:,1]

#Theoretical calculations of Mach number and pressure ratio of part 1.1
nb_theoreticalM1 = np.array([])
nb_theoreticalP1 = np.array([])
#Before the 1st throat, subsonic
for A in nb_AAt[np.where(nb_mm<x_1st_throat)]:
    out = flowtools.flowisentropic2(gamma,A,'sub')
    nb_theoreticalM1 = np.append(nb_theoreticalM1, out[0])
    nb_theoreticalP1 = np.append(nb_theoreticalP1, out[2])
#After the 1st throat, supersonic
for A in nb_AAt[np.where(nb_mm>x_1st_throat)]:
    out = flowtools.flowisentropic2(gamma,A,'sup')
    nb_theoreticalM1 = np.append(nb_theoreticalM1, out[0])
    nb_theoreticalP1 = np.append(nb_theoreticalP1, out[2])

#Calculate Mach number based on the measured pressure ratio of part 1.1
nb_calculatedM1 = np.array([])
for P in tf1_ppt:
    out = flowtools.flowisentropic2(gamma,P,'pres')
    nb_calculatedM1 = np.append(nb_calculatedM1, out[0])

#Theoretical calculations of Mach number and pressure ratio of part 1.2
#Valid minimum value of ppt2_inlet is 0.775626, otherwise at throat would be sonic
ppt2_inlet = tf2_ppt[0]
AtAsonic = flowtools.flowisentropic2(gamma, ppt2_inlet, 'pres')[-1] * (1/nb_AAt[0])
#Convert the local geometric area ratio to the local sonic area ratio 
converted_nb_AAt = AtAsonic * nb_AAt
nb_theoreticalM2 = np.array([])
nb_theoreticalP2 = np.array([])
for A in converted_nb_AAt:
    out = flowtools.flowisentropic2(gamma,A,'sub')
    nb_theoreticalM2 = np.append(nb_theoreticalM2, out[0])
    nb_theoreticalP2 = np.append(nb_theoreticalP2, out[2])

#Calculate Mach number based on the measured pressure ratio of part 1.2
nb_calculatedM2 = np.array([])
for P in tf2_ppt:
    out = flowtools.flowisentropic2(gamma,P,'pres')
    nb_calculatedM2 = np.append(nb_calculatedM2, out[0])

#Cross section area ratio of the two locations of shock in the diffuser
fd_AAt_Mloc = np.array([DiffuserAreaRatio(x_shock4, M4, M5), DiffuserAreaRatio(x_shock3, M4, M5)])
#Downstream locations of shock 4 and 3
shock_loc = np.array([x_shock4, x_shock3])

#Theoretical calculations of pressure ratio of part 2.3(1), p/pt,3
ppt3_diffuser_inlet = tf3_ppt[idx_fd_begin]
AtAsonic = flowtools.flowisentropic2(gamma, ppt3_diffuser_inlet, 'pres')[-1] * (1/fd_AAt[0])
#Convert the local geometric area ratio to the local sonic area ratio 
converted_fd_AAt_Mloc = AtAsonic * fd_AAt_Mloc
fd_theoreticalP3_Mloc = np.array([flowtools.flowisentropic2(gamma,converted_fd_AAt_Mloc[0],'sub')[2],
                                  flowtools.flowisentropic2(gamma,converted_fd_AAt_Mloc[1],'sub')[2]])

#Theoretical calculations of pressure ratio of part 2.3(2), p/pt,6
fd_theoreticalP4_Mloc = np.array([flowtools.flowisentropic2(gamma,fd_AAt_Mloc[0],'sup')[2],
                                  flowtools.flowisentropic2(gamma,fd_AAt_Mloc[1],'sup')[2]])

#Theoretical calculations of pressure ratio of part 2.3(3) after the normal shock, p/pt,5
fd_theoreticalM5_Mloc = np.array([flowtools.flowisentropic2(gamma,fd_AAt_Mloc[0],'sup')[0],
                                  flowtools.flowisentropic2(gamma,fd_AAt_Mloc[1],'sup')[0]])
out1 = flowtools.flownormalshock2(gamma,fd_theoreticalM5_Mloc[0],'mach')
out2 = flowtools.flownormalshock2(gamma,fd_theoreticalM5_Mloc[1],'mach')
#P*P1=P2/P0,2; P*P1*P0=P2/P0,1
fd_theoreticalP5_Mloc = np.array([out1[2]*out1[-1]*out1[-2],
                                  out2[2]*out2[-1]*out2[-2]])

#Plot of part 1.1
fig_1 = plt.figure(figsize=(16,8))
#Pressure ratio
P_1 = fig_1.add_subplot(211)
P_1.set_xlabel("x [mm]")
P_1.set_ylabel("p/p_t [-]")
P_1.plot(nb_mm, nb_theoreticalP1, "ko-")
P_1.plot(nb_mm, tf1_ppt, "rx-")
P_1.legend(("quasi 1D theory", "measured pressure ratio"), loc="upper right")
#Mach number
M_1 = fig_1.add_subplot(212)
M_1.set_xlabel("x [mm]")
M_1.set_ylabel("M [-]")
M_1.plot(nb_mm, nb_theoreticalM1, "ko-")
M_1.plot(nb_mm, nb_calculatedM1, "rx-")
M_1.legend(("quasi 1D theory", "computed from measured pressure ratio"), loc="lower right")

#Plot of part 1.2
fig_2 = plt.figure(figsize=(16,8))
#Pressure ratio
P_2 = fig_2.add_subplot(211)
P_2.set_xlabel("x [mm]")
P_2.set_ylabel("p/p_t [-]")
P_2.plot(nb_mm, nb_theoreticalP2, "ko-")
P_2.plot(nb_mm, tf2_ppt, "rx-")
P_2.legend(("quasi 1D theory", "measured pressure ratio"), loc="lower right")
#Mach number
M_2 = fig_2.add_subplot(212)
M_2.set_xlabel("x [mm]")
M_2.set_ylabel("M [-]")
M_2.plot(nb_mm, nb_theoreticalM2, label="quasi 1D theory", linestyle="solid", marker='o', color="black")
M_2.plot(nb_mm, nb_calculatedM2, label="computed from measured pressure ratio", linestyle="solid", marker='x', color="red")
M_2.legend(loc="upper right")

#Plot of part 2
fig_3 = plt.figure(figsize=(16,4))
P_345 = fig_3.add_subplot(111)
concat = np.concatenate((tf3_ppt, tf4_ppt, tf5_ppt, fd_theoreticalP3_Mloc, fd_theoreticalP4_Mloc, fd_theoreticalP5_Mloc))
y_max = max(concat)
y_min = min(concat)
dev = (y_max-y_min)*0.2
P_345.set_ylim(y_min-dev, y_max+dev)
P_345.set_xlabel("x [mm]")
P_345.set_ylabel("p/p_t [-]")
P_345.plot(wt_mm, tf3_ppt, label="pressure ratio of measurement 3", linestyle="solid", marker='o', color="black")
P_345.plot(wt_mm, tf4_ppt, label="pressure ratio of measurement 4", linestyle="solid", marker='o', color="blue")
P_345.plot(wt_mm, tf5_ppt, label="pressure ratio of measurement 5", linestyle="solid", marker='o', color="red")
P_345.vlines(x_shock3, y_min*0.9, y_max*1.1, color='black', linestyles='solid', linewidth=0.5)
P_345.text(x_shock3, y_min-dev/2, "shock 3", horizontalalignment='center')
P_345.vlines(x_shock4, y_min*0.9, y_max*1.1, color='black', linestyles='solid', linewidth=0.5)
P_345.text(x_shock4, y_min-dev/2, "shock 4", horizontalalignment='center')
P_345.vlines(x_shock5, y_min*0.9, y_max*1.1, color='black', linestyles='solid', linewidth=0.5)
P_345.text(x_shock5, y_min-dev/2, "shock 5", horizontalalignment='center')
P_345.plot(shock_loc,fd_theoreticalP3_Mloc, label="theoretical p/pt,3", linestyle="None", marker='v', color="grey" )
P_345.plot(shock_loc,fd_theoreticalP4_Mloc, label="theoretical p/pt,6", linestyle="None", marker='^', color="grey" )
P_345.plot(shock_loc,fd_theoreticalP5_Mloc, label="theoretical p/pt,5", linestyle="None", marker='<', color="grey" )
P_345.legend(loc="upper right")

#Plot of area ratio of nozzle block
fig_nozzle = plt.figure(figsize=(16,4))
ax_nozzle = fig_nozzle.add_subplot(111)
ax_nozzle.set_xlabel("x [mm]")
ax_nozzle.set_ylabel("A/A_t [-]")
ax_nozzle.plot(nb_mm, nb_AAt, "ko-")

#Plot of area ratio of flexible diffuser
fig_diffuser = plt.figure(figsize=(16,4))
ax_diffuser = fig_diffuser.add_subplot(111)
ax_diffuser.set_xlabel("x [mm]")
ax_diffuser.set_ylabel("A/A_t [-]")
ax_diffuser.plot(fd_mm, fd_AAt, "ko-")

fig_nozzle.savefig("Nozzle_area_ratio", dpi=300)
fig_diffuser.savefig("Diffuser_area_ratio", dpi=300)
fig_1.savefig("Part11", dpi=300)
fig_2.savefig("Part12", dpi=300)
fig_3.savefig("Part2", dpi=300)

print("\n*****OJBK!*****\n")
