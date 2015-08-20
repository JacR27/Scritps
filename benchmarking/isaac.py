import yaml as ya
import sys
class isaacRun(ya.YAMLObject):
    yaml_tag = u'!isaac3Run'
    def __init__(self,name,version,comand_line_options,compile_time,system,runtime,sample_info):
        self.name = name
        self.version = version
        self.comand_line_options = comand_line_options
        self.compile_time = compile_time
        self.system = system
        self.runtime = runtime
        self.sample_info = sample_info

sample1 = "lane 4"
version1 = "03.15.07.30"
compalation1 = "Roman"
compalation2 = "gcc-4.8.3"
compalation3 = "gcc-4.7.4"
compalation4 = "gcc-5.2.0"
compalation5 = "gcc-4.7.4-no-debug ./configure --build-type Releise"
compalation6= "gcc-4.7.4 CXXFLAGS=\"-O3\" CFLAGS=\"-O3\" ./configure"
compalation7= "gcc-5.2.0 CXXFLAGS=\"-O3\" CFLAGS=\"-O3\" ./configure"
system1 = "S4"
system2 = "swap off"
system3 = "transparent huge pages"
runtime1= "default"
runtime2= "zlibCloudflare"
runtime3= "zlibIntel"
comand1 = "default"
comand2 = "--buffer-bin no"
comand3 = "--pre-allocate-bins yes"
comand4 = "--pre-expected-bgzf-ratio 0.5"
comand5 = "--pre-expected-bgzf-ratio 0.35"
comand6 = "--temp-concurrent-load 2"
comand7 = "--temp-concurrent-load 3"
comand8 = "--temp-concurrent-load 4"
comand9 = "--temp-concurrent-load 5"
comand10 = "--temp-concurrent-load 6"
comand11 = "--temp-concurrent-load 7"
comand12 = "--temp-concurrent-load 8"
comand13 = "--temp-concurrent-load 40"
comand14 = "--base-calls-format fastq-gz --base-calls-directory /mnt/141128_HSX177_0082_BH00UABBXX/Fastq"
comand15 = "--base-calls-format fastq-gz --base-calls-directory /mnt/141128_HSX177_0082_BH00UABBXX/Fastq --clusters-at-a-time 1000000"
comand16 = "--base-calls-format fastq-gz --base-calls-directory /mnt/141128_HSX177_0082_BH00UABBXX/Fastq --clusters-at-a-time 2000000"
comand17 = "--base-calls-format fastq-gz --base-calls-directory /mnt/141128_HSX177_0082_BH00UABBXX/Fastq --clusters-at-a-time 4000000"
comand18 = "--base-calls-format fastq-gz --base-calls-directory /mnt/141128_HSX177_0082_BH00UABBXX/Fastq --clusters-at-a-time 8000000"
comand19 = "--base-calls-format fastq-gz --base-calls-directory /mnt/141128_HSX177_0082_BH00UABBXX/Fastq --clusters-at-a-time 16000000"


I01 = isaacRun("I01",version1,comand1,compalation1,system1,runtime1,sample1)
I02 = isaacRun("I02",version1,comand2,compalation1,system1,runtime1,sample1)
I03 = isaacRun("I03",version1,comand3,compalation1,system1,runtime1,sample1)
I04 = isaacRun("I04",version1,comand4,compalation1,system1,runtime1,sample1)
I05 = isaacRun("I05",version1,comand5,compalation1,system1,runtime1,sample1)
I06 = isaacRun("I06",version1,comand6,compalation1,system1,runtime1,sample1)
I07 = isaacRun("I07",version1,comand7,compalation1,system1,runtime1,sample1)
I08 = isaacRun("I08",version1,comand8,compalation1,system1,runtime1,sample1)
I09 = isaacRun("I09",version1,comand9,compalation1,system1,runtime1,sample1)
I10 = isaacRun("I10",version1,comand10,compalation1,system1,runtime1,sample1)
I11 = isaacRun("I11",version1,comand11,compalation1,system1,runtime1,sample1)
I12 = isaacRun("I12",version1,comand12,compalation1,system1,runtime1,sample1)
I13 = isaacRun("I13",version1,comand1,compalation2,system1,runtime1,sample1)
I14 = isaacRun("I14",version1,comand1,compalation3,system1,runtime1,sample1)
I15 = isaacRun("I15",version1,comand1,compalation4,system1,runtime1,sample1)
I16 = isaacRun("I16",version1,comand1,compalation1,system1,runtime2,sample1)
I17 = isaacRun("I17",version1,comand1,compalation1,system1,runtime3,sample1)
I18 = isaacRun("I18",version1,comand1,compalation5,system1,runtime1,sample1)
I19 = isaacRun("I19",version1,comand2,compalation3,system1,runtime1,sample1)
I20 = isaacRun("I20",version1,comand3,compalation3,system1,runtime1,sample1)
I21 = isaacRun("I21",version1,comand4,compalation3,system1,runtime1,sample1)
I22 = isaacRun("I22",version1,comand5,compalation3,system1,runtime1,sample1)
I23 = isaacRun("I23",version1,comand6,compalation3,system1,runtime1,sample1)
I24 = isaacRun("I24",version1,comand7,compalation3,system1,runtime1,sample1)
I25 = isaacRun("I25",version1,comand8,compalation3,system1,runtime1,sample1)
I26 = isaacRun("I26",version1,comand13,compalation3,system1,runtime1,sample1)
I27 = isaacRun("I27",version1,comand1,compalation6,system1,runtime1,sample1)
I28 = isaacRun("I28",version1,comand1,compalation1,system2,runtime1,sample1)
I29 = isaacRun("I29",version1,comand1,compalation1,system3,runtime1,sample1)
I30 = isaacRun("I30",version1,comand1,compalation3,system1,runtime2,sample1)
I31 = isaacRun("I31",version1,comand1,compalation3,system1,runtime3,sample1)
I32 = isaacRun("I32",version1,comand14,compalation1,system1,runtime1,sample1)
I33 = isaacRun("I33",version1,comand1,compalation7,system1,runtime1,sample1)
I34 = isaacRun("I34",version1,comand1,compalation3,system3,runtime1,sample1)
I35 = isaacRun("I35",version1,comand1,compalation3,system2,runtime1,sample1)
I36 = isaacRun("I36",version1,comand15,compalation1,system1,runtime1,sample1)
I37 = isaacRun("I37",version1,comand16,compalation1,system1,runtime1,sample1)
I38 = isaacRun("I38",version1,comand17,compalation1,system1,runtime1,sample1)
I39 = isaacRun("I39",version1,comand18,compalation1,system1,runtime1,sample1)
I40 = isaacRun("I40",version1,comand19,compalation1,system1,runtime1,sample1)


print(ya.dump(eval(sys.argv[1]),default_flow_style=False))
