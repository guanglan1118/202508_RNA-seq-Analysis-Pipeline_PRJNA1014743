bsub -P job_name \
     -J job_name \
     -n 2 \
     -R "rusage[mem=8GB]" \
     -eo job_name.err \
     -oo job_name.out \
     "sh TS_linker_workflow.sh"


#!/bin/bash
#BSUB -P job_name            # 项目名（集群用来记账或分配资源）
#BSUB -J job_name            # 任务名（方便用 bjobs 查看）
#BSUB -n 2                   # 使用 2 个 CPU 核
#BSUB -R "rusage[mem=8GB]"   # 每个核 8GB 内存
#BSUB -eo job_name.err       # 错误日志输出到 job_name.err
#BSUB -oo job_name.out       # 标准输出到 job_name.out
#BSUB -W 04:00               # (可选) 最大运行时间 4 小时

echo "Job started on $(date)"
echo "Running on $(hostname)"

# 执行你的工作流程脚本
sh TS_linker_workflow.sh

echo "Job finished on $(date)"

~~~
#!/bin/bash
#BSUB -P salmon_index_build        # 项目名（按需改成你们集群要求的）
#BSUB -J salmon_index              # 任务名
#BSUB -n 8                         # 请求 8 个 CPU 核
#BSUB -R "rusage[mem=8GB]"         # 每个核 8GB 内存 (总共 64GB)
#BSUB -W 12:00                     # 最大运行时间 12 小时 (按需要调整)
#BSUB -oo salmon_index.out         # 标准输出日志
#BSUB -eo salmon_index.err         # 错误日志

echo "Job started on $(date)"
echo "Running on node $(hostname)"

# 切换到工作目录 (替换为你放参考文件的路径)
cd /research/groups/yanggrp/home/glin/work_2025/Sep/project_PRJNA1014743/ref

# 运行 salmon index
salmon index \
  -t gencode.v44.gentrome.fa \
  -d decoys.txt \
  -i salmon_gencode_v44_decoy \
  --gencode \
  -p 8

echo "Job finished on $(date)"

~~~
在命令行运行：
bsub < salmon_index.lsf

[info] Building perfect hash
[info] Index built successfully





