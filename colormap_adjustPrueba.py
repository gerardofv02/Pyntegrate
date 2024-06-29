import Pyntegrate
from matplotlib import pyplot as plt

fig = plt.figure()
# fig = Pyntegrate.colormap_adjust.color_test("red")
# plt.show()

fig = Pyntegrate.colormap_adjust.smart_colormap(0,100)
plt.show()