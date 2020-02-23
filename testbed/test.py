import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
def exponential_decay(t, y): return -0.5 * y
sol = solve_ivp(exponential_decay, [0, 10], [2, 4, 8]);
plt.figure()
plt.plot(sol.t, sol.y[0])
plt.savefig("test.png")
