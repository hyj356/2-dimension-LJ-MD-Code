# 编译方法

​		确保你的电脑上面安装了make和gfortran, 然后进入src目录下面, 输入:

```shell
make
```

​		那么默认就会开启openMP并行进行编译. 并获得一个名为target的可执行文件, 如果不想开启openMP并行, 可以输入:

```shell
make serial
```

​		那么就会获得一个名为target_serial的可执行文件.

​		通过输入:

```shell
export OMP_NUM_THREADS=5
```

​		可以将openMP的并行核数设置为5核, 请同学们根据自己的硬件情况进行设置.



# 使用方法

​		获得了可执行文件之后, 在src目录下直接输入可执行文件的名字, 程序就会自动读取parameter.nml文件中的参数并以文件中的参数进行MD仿真. 这里说明一下parameter.nml文件的名字是不能改的, 只有里面的内容可以修改. parameter.nml文件的内容和物理意义如下:

```
&InputVar
deltaT      =      0.002d0        ! 时间步长
density     =      0.8d0          ! 原子分布的密度
nx          =      50             ! x方向上晶格个数
ny          =      50             ! y方向上晶格个数
stepAvg     =      100            ! 每隔多少步平均一次物理量
stepEquil   =      0              ! 目前仿真中暂时没用到的物理量
stepLimit   =      5000           ! 一共循环多少步
temperature =      1.d0           ! 模拟在什么温度下进行
/

```

​		请同学们根据自己的需要设置对应的仿真参数.