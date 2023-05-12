# 先在R里面把exp 导出
```R
exp = as.matrix(as.data.frame(scRNA@assays$RNA@data))
write.csv(exp, "xx.csv") # 把counts导出
```
# 再在python里面获取
```python
import pandas as pd
from scipy.stats import spearmanr
df = pd.read_csv("xx.csv").T # 读取counts
spearmanr(df.loc[:, 'gene1'], df.loc[:, 'gene2']) # 得到相关性
```