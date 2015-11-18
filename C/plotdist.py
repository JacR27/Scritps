import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def plotdist():
    data = [45.9,0.3,5.6,2.6,3,6,9.2,27.5]
    pos = [i for i in range(len(data))]
    scores = [0,7,11,22,27,32,37,42]
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    barWidth = .75
    ax1.bar(pos,data,width=barWidth)
    ax1.set_xticks([(i+(barWidth/2)) for i in pos])
    ax1.set_xticklabels(scores)
    ax1.set_ylabel("Percentage of base-calls")
    ax1.set_xlabel("Quality score")
    ax1.set_title("Q score distribution ~30x genome") 
    
    fig.savefig("dist.png")




plotdist()
