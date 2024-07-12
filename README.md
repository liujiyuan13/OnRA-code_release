# On the Correlation Measurement of Data Representations (数据表示的相关性度量方法)

Matlab implementation for CJC (计算机学报) paper:

- Jiyuan Liu, Xinwang Liu, Zhiping Cai, En Zhu, Rusheng Ju: [On the Correlation Measurement of Data Representations (数据表示的相关性度量方法)](http://cjc.ict.ac.cn/online/onlinepaper/ljy-202479162037.pdf). Chinese Journal of Computers (计算机学报), CJC, 2024.

## Introduction
**Abstract**

The correlation measurement between data representations is the basis of machine learning and artificial intelligence techniques, However, existing methods either overlook the global information on the involved representations partially or scale poorly, To address the issue, this paper proposes the Representation Alignment (RA) to globally measure the correlation between two arbitrary data representations in linear time to the dimension and the number of samples, We further extend the representation alignment to the Contrastive Representation Alignment (CRA), where the concerned representations are two augmentations (in the scope of contrastive learning) of data samples. In addition, the two alignments can be naturally applied in multi-view learning systems by maximizing them on the latent representations of different data views so as to integrate their complementary information, To validate this, we develop two novel multi-view clustering algorithms and achieve state-of-the-art performance on seven benchmark datasets.

数据表示之间的相关性度量是机器学习和人工智能技术的基石。
然而现有的度量方法要么数据表示的全局信息考虑不足，要么复杂度较高，限制了相关技术的进一步发展。
为解决上述问题，本文提出一种数据表示的对齐度量方法，称为表示对齐(Representation Alignment, RA)。
此度量方法能够全局性地衡量任意两个数据表示之间的相关性，且其在样本数量和特征维度上的计算复杂度均为线性。
在此基础上，我们将RA扩展到了对比学习领域，进一步提出了基于对比的表示对齐(Contrastive Representation Alignment, CRA)度量方法。
上述两个度量方法可自然地用于多视图学习场景，即可通过最大化不同视图数据之间的RA和CRA来融合各个视图之间的信息。
为验证这一点，我们还提出了两个新颖的多视图聚类算法，并在七个基准数据集上取得了领先的聚类性能。

## Code structure

```
...
+ eval				# Matlab functions for evaluation
+ plot 				# plot results
.gitignore 			
EuDist2.m			# tool function
get_res.m    		# obtain results
LICENSE.py 			# license file
mc_cra.m            # MCCRA algorithm
mc_ra.m             # MCRA algorithm
README.md 
run_file_cra.m 	    # run file (example) of MCCRA algorithm
run_file_ra.m 		# run file (example) of MCRA algorithm
```

## Usage

1. Clone to the local.
```
> git clone https://github.com/liujiyuan13/OnRA-code_release.git OnRA-code_release
```
2. Run the algorithms.
```
> run_file_ra
> run_file_cra
```
3. Get results. 
```
> get_res			# 'res_out' format: [dataset, metrics(acc,nmi,purity), algorithm(MCKA, MCRA, MCCRA)]
```


## Citation

If you find our code useful, please cite:

    @article{liu2024onrachi,
        title        = {数据表示的相关性度量方法},
        author       = {刘吉元 and 刘新旺 and 蔡志平 and 祝恩 and 鞠儒生},
        year         = 2024,
        journal      = {计算机学报},
        volume       = {第47卷},
        pages        = {1568--1581},
        issue        = {第7期}
    }
    @article{liu2024onraen,
        title        = {On the Correlation Measurement of Data Representations},
        author       = {Jiyuan Liu and Xinwang Liu and Zhiping Cai and En Zhu and Rusheng Ju},
        year         = 2024,
        journal      = {Chinese Journal of Computers},
        volume       = 47,
        pages        = {1568--1581},
        issue        = 7
    }


## Licence

This repository is under [GPL V3](https://github.com/liujiyuan13/OnRA-code_release/blob/main/LICENSE).

## More
- For more related researches, please visit my homepage: [https://liujiyuan13.github.io](https://liujiyuan13.github.io).
- For data and discussion, please message me: liujiyuan13@nudt.edu.cn




