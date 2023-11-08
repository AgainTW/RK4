#不同能量對繩波影響，還沒做
import math
import numpy as np
import os
from scipy.special import gamma
from matplotlib import pyplot as plt

def PDE(y,coefficient_1,coefficient_2,coefficient_3,coefficient_4,n):
	for j in range(1):
		for i in range(1,100,1):
			y[n+1][i]=coefficient_1*y[n][i]+coefficient_2*(y[n][i+1]+y[n][i-1]-2*y[n][i])-coefficient_3*y[n-1][i]+coefficient_4
	return y

#數值設定
###一般波函數部分
step=120000
y_1=np.zeros(101*step)
y_2=np.zeros(101*step)
y_3=np.zeros(101*step)
y_4=np.zeros(101*step)
y_1=y_1.reshape((step,101))
y_2=y_2.reshape((step,101))
y_3=y_3.reshape((step,101))
y_4=y_4.reshape((step,101))
y_tmp_1=np.zeros(step)
y_tmp_2=np.zeros(step)
y_tmp_3=np.zeros(step)
y_tmp_4=np.zeros(step)
tmp_1=[]
tmp_2=[]
tmp_3=[]
tmp_4=[]
period_1=[]
period_2=[]
period_3=[]
period_4=[]

l=1						# 弦長
Tau=40					# 張力常數
p=0.01					# 線密度
v=math.sqrt(Tau/p)		# 波速
t=12.0					# 經過時間,三次回到原位
b_1=1					# y1摩擦力項常數
b_2=1					# y2摩擦力項常數
b_3=1					# y3摩擦力項常數
b_4=1					# y4摩擦力項常數
c_1=0					# y1外力項常數
c_2=0					# y2外力項常數
c_3=0					# y3外力項常數
c_4=0					# y4摩擦力項常數
zto=np.linspace(0, 1, 101)
ztt=np.linspace(0, 12, step)
###複合變數部分
delta_t=t/step 			# delta_t=0.0001時，較為精確
delta_x=l/99
#沒摩擦力沒重力的
coeffi_1=b_1*delta_t+1
coefficient_1_1=(b_1*delta_t+2)/coeffi_1
coefficient_1_2=(v*v*delta_t*delta_t)/(delta_x*delta_x*coeffi_1)
coefficient_1_3=1/coeffi_1
coefficient_1_4=c_1*delta_t*delta_t/coeffi_1
#有摩擦力沒重力的
coeffi_2=b_2*delta_t+1
coefficient_2_1=(b_2*delta_t+2)/coeffi_2
coefficient_2_2=(v*v*delta_t*delta_t)/(delta_x*delta_x*coeffi_2)
coefficient_2_3=1/coeffi_2
coefficient_2_4=c_2*delta_t*delta_t/coeffi_2
#有摩擦力有重力的
coeffi_3=b_3*delta_t+1
coefficient_3_1=(b_3*delta_t+2)/coeffi_3
coefficient_3_2=(v*v*delta_t*delta_t)/(delta_x*delta_x*coeffi_3)
coefficient_3_3=1/coeffi_3
coefficient_3_4=c_3*delta_t*delta_t/coeffi_3
#沒摩擦力有重力的
coeffi_4=b_4*delta_t+1
coefficient_4_1=(b_4*delta_t+2)/coeffi_4
coefficient_4_2=(v*v*delta_t*delta_t)/(delta_x*delta_x*coeffi_4)
coefficient_4_3=1/coeffi_4
coefficient_4_4=c_4*delta_t*delta_t/coeffi_4

print('set')
#初始化
###y_1:
for i in range(1,100,1):
	if(i<=80):
		y_1[0][i]=0.125*i/100
	else:
		y_1[0][i]=0.5*(100-i)/100
y_1[0][100]=0
for i in range(0,101,1):
	y_1[1][i]=y_1[0][i]
###y_2:
for i in range(1,100,1):
	if(i<=60):
		y_2[0][i]=0.125*i/100
	else:
		y_2[0][i]=0.5*(100-i)/100
y_2[0][100]=0
for i in range(0,101,1):
	y_2[1][i]=y_2[0][i]
###y_3:
for i in range(1,100,1):
	if(i<=40):
		y_3[0][i]=0.125*i/100
	else:
		y_3[0][i]=0.5*(100-i)/100
y_3[0][100]=0
for i in range(0,101,1):
	y_3[1][i]=y_3[0][i]
###y_4:
for i in range(1,100,1):
	if(i<=20):
		y_4[0][i]=0.125*i/100
	else:
		y_4[0][i]=0.5*(100-i)/100
y_4[0][100]=0
for i in range(0,101,1):
	y_4[1][i]=y_4[0][i]

print('initialization')

#計算
for i in range(1,step-1,1):
	y_1=PDE(y_1,coefficient_1_1,coefficient_1_2,coefficient_1_3,coefficient_1_4,i)
print('Finsh 1')
for i in range(1,step-1,1):	
	y_2=PDE(y_2,coefficient_2_1,coefficient_2_2,coefficient_2_3,coefficient_2_4,i)
print('Finsh 2')
for i in range(1,step-1,1):
	y_3=PDE(y_3,coefficient_3_1,coefficient_3_2,coefficient_3_3,coefficient_3_4,i)
print('Finsh 3')
for i in range(1,step-1,1):
	y_4=PDE(y_4,coefficient_4_1,coefficient_4_2,coefficient_4_3,coefficient_4_4,i)
print('Finsh 4')

#找波峰
for i in range(1,step-1,1):
	y_tmp_1[i-1]=y_1[i][75]
for i in range(1,step-1,1):
	y_tmp_2[i-1]=y_2[i][75]
for i in range(1,step-1,1):
	y_tmp_3[i-1]=y_3[i][75]
for i in range(1,step-1,1):
	y_tmp_4[i-1]=y_4[i][75]

#只取極值
judge=0
tmp_1.append(y_tmp_1[0])
for i in range(1,step-1,1):
	if(y_tmp_1[i]<y_tmp_1[i+1] and judge==0):
		tmp_1.append(y_tmp_1[i])
		judge=1
	elif(y_tmp_1[i]>y_tmp_1[i+1] and judge==1):
		tmp_1.append(y_tmp_1[i])
		period_1.append(i/10000)
		judge=0
judge=0
tmp_2.append(y_tmp_2[0])
for i in range(1,step-1,1):
	if(y_tmp_2[i]<y_tmp_2[i+1] and judge==0):
		tmp_2.append(y_tmp_2[i])
		judge=1
	elif(y_tmp_2[i]>y_tmp_2[i+1] and judge==1):
		tmp_2.append(y_tmp_2[i])
		period_2.append(i/10000)
		judge=0
judge=0
tmp_3.append(y_tmp_3[0])
for i in range(1,step-1,1):
	if(y_tmp_3[i]<y_tmp_3[i+1] and judge==0):
		tmp_3.append(y_tmp_3[i])
		judge=1
	elif(y_tmp_3[i]>y_tmp_3[i+1] and judge==1):
		tmp_3.append(y_tmp_3[i])
		period_3.append(i/10000)
		judge=0
judge=0
tmp_4.append(y_tmp_4[0])
for i in range(1,step-1,1):
	if(y_tmp_4[i]<y_tmp_4[i+1] and judge==0):
		tmp_4.append(y_tmp_4[i])
		judge=1
	elif(y_tmp_4[i]>y_tmp_4[i+1] and judge==1):
		tmp_4.append(y_tmp_4[i])
		period_4.append(i/10000)
		judge=0

#週期計算
for i in range(len(period_1)-1,0,-1):
	period_1[i]=period_1[i]-period_1[i-1]
for i in range(len(period_2)-1,0,-1):
	period_2[i]=period_2[i]-period_2[i-1] 
for i in range(len(period_3)-1,0,-1):
	period_3[i]=period_3[i]-period_3[i-1] 
for i in range(len(period_4)-1,0,-1):
	period_4[i]=period_4[i]-period_4[i-1] 

ztp_1=np.linspace(0, 12, len(period_1))
ztp_2=np.linspace(0, 12, len(period_2))
ztp_3=np.linspace(0, 12, len(period_3))
ztp_4=np.linspace(0, 12, len(period_4))
zttmp_1=np.linspace(0, 12, len(tmp_1))
zttmp_2=np.linspace(0, 12, len(tmp_2))
zttmp_3=np.linspace(0, 12, len(tmp_3))
zttmp_4=np.linspace(0, 12, len(tmp_4))

#固定某點作圖
print('start plot')
for i in range(1):
	plt.figure()
	plt.title("step=10000,t=1.0") 		
	plt.xlim((0, 12))
#	plt.ylim((-10.0, 10.0))
	l1, =plt.plot(zttmp_1, tmp_1,'r.',label='x=80')
	l2, =plt.plot(zttmp_2, tmp_2,'y.',label='x=60')
	l3, =plt.plot(zttmp_3, tmp_3,'g.',label='x=40')
	l4, =plt.plot(zttmp_4, tmp_4,'b.',label='x=20')
	plt.legend(loc='upper right')
	plt.savefig('C:/Users/阿甘/Desktop/祭物3/photo_24/不同波峰位置作圖.png')
	plt.close()

#週期作圖
for i in range(1):
	plt.figure()
	plt.title("step=10000,t=1.0") 		
	plt.xlim((0, 12))
	l1, =plt.plot(ztp_1, period_1,'r.',label='x=80')
	plt.legend(loc='upper right')
	plt.savefig('C:/Users/阿甘/Desktop/祭物3/photo_24/x=80.png')
	plt.close()
for i in range(1):
	plt.figure()
	plt.title("step=10000,t=1.0") 		
	plt.xlim((0, 12))
	l2, =plt.plot(ztp_2, period_2,'y.',label='x=60')
	plt.legend(loc='upper right')
	plt.savefig('C:/Users/阿甘/Desktop/祭物3/photo_24/x=60.png')
	plt.close()
for i in range(1):
	plt.figure()
	plt.title("step=10000,t=1.0") 		
	plt.xlim((0, 12))
	l2, =plt.plot(ztp_3, period_3,'g.',label='x=40')
	plt.legend(loc='upper right')
	plt.savefig('C:/Users/阿甘/Desktop/祭物3/photo_24/x=40.png')
	plt.close()
for i in range(1):
	plt.figure()
	plt.title("step=10000,t=1.0") 		
	plt.xlim((0, 12))
	l2, =plt.plot(ztp_4, period_4,'b.',label='x=80')
	plt.legend(loc='upper right')
	plt.savefig('C:/Users/阿甘/Desktop/祭物3/photo_24/x=20.png')
	plt.close()

'''
#動畫作圖
print('start plot')
for i in range(1,step,83):
	plt.figure()
	plt.title("step=120000,t=12.0,b=1,c=0") 		
	plt.xlim((0, 1))
	plt.ylim((-10.0, 10.0))
	l1, =plt.plot(zto, y_1[i],'r',label='max_amp=0.1')
	l2, =plt.plot(zto, y_2[i],'yv',label='max_amp=1')
	l3, =plt.plot(zto, y_3[i],'g^',label='max_amp=5')
	l4, =plt.plot(zto, y_4[i],'b.',label='max_amp=10')
	plt.legend(loc='upper right')
	plt.savefig('C:/Users/阿甘/Desktop/祭物3/photo_22/'+str(i)+'.png')
	plt.close()
'''

