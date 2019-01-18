def print_hamming(coordinates):
	A,B,C,D,E,F= coordinates
	width = max(2*max(A,C)+F,2*max(B,D)+F) + 16
	#Print top portion of plot
	x=min(A,C); y=max(A,C)
	while (y>x and y>1):
		if A>C:
			print(('\\'+ (2*(y-1)+F)*(' ')+ ' ').center(width))
		if C>A:
			print((' '+ (2*(y-1)+F)*(' ')+ '/').center(width))
		y-=1
	while (y==x and y>1):
		print(('\\'+ (2*(y-1)+F)*(' ')+ '/').center(width))
		y-=1; x-=1
	if y==1:
		if x==1:
			print(('\\'+ (2*(y-1)+F)*('_')+ '/').center(width)); y-=1; x-=1
		elif C>A:
			print((' '+ (2*(y-1)+F)*('_')+ '/').center(width)); y-=1
		elif A>C:
			print(('\\'+ (2*(y-1)+F)*('_')+ ' ').center(width)); y-=1
	else:
		print((' '+ (2*(y-1)+F)*('_')+ ' ').center(width))
	#Print middle portion of plot
	if E > 0:
		print((('|'+F*' '+'|').center(width)+'\n')*(E-1)  +  ('|'+F*'_'+'|').center(width))
			
	#Print bottom portion of plot
	X=min(B,D); Y=max(B,D)
	while x<X:
		print(('/'+ (2*y+F)*(' ')+ '\\').center(width))
		x+=1; y+=1
	while (x==X and y<Y):
		if B>D:
			print(('/'+ (2*y+F)*(' ')+ ' ').center(width)); y+=1
		if D>B:
			print((' '+ (2*y+F)*(' ')+ '\\').center(width)); y+=1
	print('')
##Horizontal Extension Printing##
def print_hamming_horizontally(coordinates):
	A,B,C,D,E,F= coordinates
	print((A*'_'+' '+F*'_'+' '+C*'_').rjust(C+F+2+B)+'\n'\
		+(('|'+F*' '+'|').rjust(F+2+max(A,B))+'\n')*(E-1)\
		+(B*'_'+'|'+F*'_'+'|'+D*'_').rjust(D+F+2+A)+'\n')

