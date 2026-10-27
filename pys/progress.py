# n = 789

# prev = 0
# for i in range(n):
# 	frac = int(100*i/n)
# 	if frac > 0 and frac%10 == 0 and frac > prev:
# 		print(frac,'%')#,i/n,round(i/n,2))
# 		prev = frac
# 		# print(str(100*i/n),"percent done")

# 	# else:
# 	# 	print('===')

# 	if i == n-1:
# 		print('100%')



def progress(i, n, prev):
	frac = int(10000*i/n)
	if frac > 0 and frac > prev and frac%100 == 0:
		return frac/100