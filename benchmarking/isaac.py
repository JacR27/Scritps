import yaml as ya

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
system1 = "S4"
runtime1= "default"
comand1 = "default"
comand2 = "--buffer-bin no"
comand3 = "--pre-allocate-bins yes"

I1 = isaacRun("I1",version1,comand1,compalation1,system1,runtime1,sample1)
I2 = isaacRun("I2",version1,comand2,compalation1,system1,runtime1,sample1)
I3 = isaacRun("I3",version1,comand3,compalation1,system1,runtime1,sample1)
print(ya.dump(I1,default_flow_style=False))
