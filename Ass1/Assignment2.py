import requests
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from bs4 import BeautifulSoup

data = np.genfromtxt('hygdata_v3.csv', delimiter=',',dtype='str')
ra = data[1:-1,23]
ra = ra.astype(np.float)
dec = data[1:-1,24]
dec = dec.astype(np.float)
mag = data[1:-1,13]
mag = mag.astype(np.float)
typ = data[1:-1,15]
var = data[1:-1,34]
ind1 = np.where(var=='')
ind2 = np.where(var!='')

F = 100*np.exp(mag/(-2.5))

ra_O = []
ra_B = []
ra_F = []
ra_M = []

dec_O = []
dec_B = []
dec_F = []
dec_M = []

mag_O = []
mag_B = []
mag_F = []
mag_M = []

for i in range(len(typ)):
    try:
        if(typ[i][0]=='O'):
            ra_O.append(ra[i])
            dec_O.append(dec[i])
            mag_O.append(mag[i])
        elif(typ[i][0]=='B'):
            ra_B.append(ra[i])
            dec_B.append(dec[i])
            mag_B.append(mag[i])
        elif(typ[i][0]=='F'):
            ra_F.append(ra[i])
            dec_F.append(dec[i])
            mag_F.append(mag[i])
        elif(typ[i][0]=='M'):
            ra_M.append(ra[i])
            dec_M.append(dec[i])
            mag_M.append(mag[i])
    except:
        continue

F_O = 500*np.exp(np.array([mag_O])/(-2.5))
F_B = 500*np.exp(np.array([mag_B])/(-2.5))
F_F = 500*np.exp(np.array([mag_F])/(-2.5))
F_M = 500*np.exp(np.array([mag_M])/(-2.5))

ra_O = np.array([ra_O]) - np.pi
ra_B = np.array([ra_B]) - np.pi
ra_F = np.array([ra_F]) - np.pi
ra_M = np.array([ra_M]) - np.pi

fig = plt.figure(figsize=(20, 10),constrained_layout=True)

gs = fig.add_gridspec(2, 2)

axO = fig.add_subplot(gs[0,0], projection="mollweide")
axO.set_facecolor('black')
axO.set_title('O',color = 'black',fontsize='xx-large')
axO.scatter(ra_O,dec_O,s=F_O,c='white')
axO.xaxis.set_ticks_position('top')
tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
tick_labels = np.remainder(tick_labels+360+180,360)
axO.set_xticklabels(tick_labels) 
for tick in axO.xaxis.get_ticklabels():
    tick.set_fontsize('large')
    tick.set_color('white')
    tick.set_weight('bold')
    
axB = fig.add_subplot(gs[0,1], projection="mollweide")
axB.set_facecolor('black')
axB.set_title('B',color = 'black',fontsize='xx-large')
axB.scatter(ra_B,dec_B,s=F_B/10,c='white')
axB.xaxis.set_ticks_position('top')
tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
tick_labels = np.remainder(tick_labels+360+180,360)
axB.set_xticklabels(tick_labels) 
for tick in axB.xaxis.get_ticklabels():
    tick.set_fontsize('large')
    tick.set_color('white')
    tick.set_weight('bold')

axF = fig.add_subplot(gs[1,0], projection="mollweide")
axF.set_facecolor('black')
axF.set_title('F',color = 'black',fontsize='xx-large')
axF.scatter(ra_F,dec_F,s=F_F/50,c='white')
axF.xaxis.set_ticks_position('top')
tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
tick_labels = np.remainder(tick_labels+360+180,360)
axF.set_xticklabels(tick_labels) 
for tick in axF.xaxis.get_ticklabels():
    tick.set_fontsize('large')
    tick.set_color('white')
    tick.set_weight('bold')
    
axM = fig.add_subplot(gs[1,1], projection="mollweide")
axM.set_facecolor('black')
axM.set_title('M',color = 'black',fontsize='xx-large')
axM.scatter(ra_M,dec_M,s=F_M/5,c='white')
axM.xaxis.set_ticks_position('top')
tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
tick_labels = np.remainder(tick_labels+360+180,360)
axM.set_xticklabels(tick_labels) 
for tick in axM.xaxis.get_ticklabels():
    tick.set_fontsize('large')
    tick.set_color('white')
    tick.set_weight('bold')

plt.show()

page = requests.get("http://astrosat.iucaa.in/czti/?q=grb")
soup = BeautifulSoup(page.content, 'html.parser')

data2 = []
for i in soup.find_all('tr')[1:-1]:
    try:
        temp = i.find_all('td')[3].get_text().strip('\n\t\t\t\t').split(',')
        temp = np.array(temp)
        temp = temp.astype(np.float)
    except:
        continue
    data2.append(temp)

data2 = np.array(data2)
data2 = data2*np.pi/180
data2[:,0] = data2[:,0] - np.pi
ra = np.array(ra) - np.pi

fig2 = plt.figure(figsize=(20, 10),constrained_layout=True)

ax = fig2.add_subplot(111, projection="mollweide")
ax.set_facecolor('black')
ax.set_title('GRB',color = 'black',fontsize='xx-large')
ax.scatter(ra[ind1[0]],dec[ind1[0]],s=F/100,c='white')
ax.scatter(data2[:,0],data2[:,1],s=10,c='blue')
ax.scatter(ra[ind2[0]],dec[ind2[0]],s=10,c='yellow')
tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
tick_labels = np.remainder(tick_labels+360+180,360)
ax.set_xticklabels(tick_labels) 
for tick in ax.xaxis.get_ticklabels():
    tick.set_fontsize('large')
    tick.set_color('white')
    tick.set_weight('bold')
plt.show()
