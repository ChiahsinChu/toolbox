# make plots

## set style

```python
from toolbox import plot
import numpy as np
import matplotlib.pyplot as plt

plt.plot(np.arange(10), np.arange(10), color="blue")
plot.use_style("pub")
plt.plot(np.arange(10), np.arange(10), color="blue")

plt.savefig("mpl_style_example.png")

plt.show()
```

Output:
![mpl_style_example](figures/mpl_style_example.png)
