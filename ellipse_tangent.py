from sympy import *
from math import sin, cos, radians, sqrt
import matplotlib.pyplot as plt

A,B,X,Y,x,y,a,b = symbols('A,B,X,Y,x,y,a,b')

equation = [x**2/A**2 + y**2/B**2 - 1, Y - a*X - b, y - a*x - b]
equation = [x**2/A**2 + (a*x+b)**2/B**2 - 1, Y - a*X - b]
equation = [x**2/A**2 + (a*x+Y-a*X)**2/B**2 - 1]
# discriminant must be 0 because we want the lines to intersects ellipse only once
equation.append(discriminant(equation[0], (x,a)))
equation.append(y-a*x-Y+a*X)
#for e in equation: print(e)
#for sol in solve(equation, (x,y,a), dict=True): print(sol)

A = 2
B = 3
X = 4
Y = 3

x1 = A**2*(B**2*X + Y*sqrt(-A**2*B**2 + A**2*Y**2 + B**2*X**2))/(A**2*Y**2 + B**2*X**2)
y1 = B**2*(A**2*Y - X*sqrt(-A**2*B**2 + A**2*Y**2 + B**2*X**2))/(A**2*Y**2 + B**2*X**2)

x2 = A**2*(B**2*X - Y*sqrt(-A**2*B**2 + A**2*Y**2 + B**2*X**2))/(A**2*Y**2 + B**2*X**2)
y2 = B**2*(A**2*Y + X*sqrt(-A**2*B**2 + A**2*Y**2 + B**2*X**2))/(A**2*Y**2 + B**2*X**2)

plt.plot([A*cos(radians(t)) for t in range(0, 360)], [B*sin(radians(t)) for t in range(0, 360)])
plt.plot([X-0.1,X+0.1], [Y-0.1,Y+0.1])
plt.plot([X+0.1,X-0.1], [Y-0.1,Y+0.1])
plt.plot([x1-0.1,x1+0.1], [y1-0.1,y1+0.1])
plt.plot([x1+0.1,x1-0.1], [y1-0.1,y1+0.1])
plt.plot([x2+0.1,x2-0.1], [y2-0.1,y2+0.1])
plt.plot([x2-0.1,x2+0.1], [y2-0.1,y2+0.1])
plt.show()

