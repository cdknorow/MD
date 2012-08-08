## \package MD.analysis.reconstruct_lattice 
# \brief This module is used to reconstruct the real space lattice from a
# structure factor q vectors(note this isn't a very reliable method)
# 
# 

##\ brief Find the real space basis
# 
# This method is primitive and unreliable and will most likley only work for
# simple lattice structures which have only 3 basis vectors
#
# This function also writes out qnumbers.txt which will be an organized list of
# the |q| as well recipricol.txt which holds the recipricol lattice vectors 
#
# \returns real space lattice vectors
# \returns basis vectors in k space
#
# \param QQ |q| where q is the primitive vectors in k space
# \param PM primitive vectors q in k space
def reconstruct_lattice(QQ,PM):
	#delete the 0 value from Q
	d=QQ.index(0)
	del(QQ[d])
	del(PM[d])
	#Find the indexes of the min values in QQ that correspond to the first
	#peaks of the fcc crystal
	#Then add the PM index to b
	b=[]
	q=[]
	fid = open('qnumber.txt','w')
	PP=copy.deepcopy(QQ)
	PP.sort()
	for i in PP:
		fid.write(('%f \n')%(i))
	try:
		for i in range(12):
			ZZ=copy.deepcopy(QQ)
			ZZ.sort()
			index = QQ.index(min(ZZ))
			b.append(PM[index])
			q.append(QQ[index])
			del(QQ[index])
			del(PM[index])
	except:
		pass
	#Find the R-space vectors
	#Organize into three rectangles of the FCC
	for i in b:
		print i,"\n"
	b1=b[0]
	for i in b:
		if i!=b1:
			b2=i
			break
	for i in b:
		if i!=b1 and i!=b2:
			b3=i
			break
	print "Found b1,b2,b3"
	print b1,"\n", b2,"\n",b3
	fid =open('recipricol.txt','w')
	fid.write(('%f %f %f \n %f %f %f \n %f %f %f\n')%(b1[0],
		b1[1],b1[2],b2[0],b2[1],b2[2],b3[0],b3[1],b3[2]))
	fid.close()
	den=np.dot(b1,np.cross(b2,b3))
	a1 = np.cross(b2,b3)
	a2 = np.cross(b3,b1)
	a3 = np.cross(b1,b2)
	a1 /= den
	a2 /= den
	a3 /= den
	a1 /= 1/(math.pi*2)
	a2 /= 1/(math.pi*2)
	a3 /= 1/(math.pi*2)

