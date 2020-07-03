##1.ggplot(), 括号内data=xxx, mapping=aes(x= ,y= )
##2.geom_xxx(), xxx包括point，bar，boxplot，density，histogram一般括号内有x,y,colour,position,shape等等
##3.scale_xxx(),常用有scale_y_log10(),scale_colour_manual(values=rainbow(7))
##4.stat_xxx(),统计变换，例如stat_smooth(),平滑曲线化
##5.coord_xxx(),坐标变换，例如coord_flip(),横纵坐标调换，或者coord_polar，极坐标化（变成pie chart）
##6.+，图层概念
##7.facet_wrap(~xxx),根据变量xxx分面
##8。theme(),ggtitle("xxx"),xlab("xxx"),ylab("xxx"),theme默认为灰色，可以调整theme_xxx()变成其他主题

require(ggplot2)
data(diamonds)
set.seed(42)
small <- diamonds[sample(nrow(diamonds),1000),]
head(small)
summary(small)

##################
##散点图，直观的展示了整个数据的分布
p <- ggplot(data = small, mapping = aes(x = carat,y = price))
p + geom_point()
p <- ggplot(data=small, mapping=aes(x=carat,y=price,shape=cut))
p + geom_point()
p <- ggplot(data=small, mapping=aes(x=carat,y=price,colour=color))
p + geom_point()

###既可以像上边那样在ggplot后定义属性，例如xy轴，形状，颜色
###也可以像下面这样，根据不同集合对象geom_xxxx,在geom后添加属性

p <- ggplot(small)
p + geom_point(aes(x=carat,y=price,shape=cut,colour=color))

##通常将不同图层共用的映射，例如xy轴提供给ggplot

##############################
##直方图geom_histogram，信息量=1变量+该变量的分布统计
##适用于离散的数据，在跨度较大范围内进行分类统计
ggplot(small)+geom_histogram(aes(x=price))

##直方图也可以将1个变量进行不同的填充，比如价格对应的切工

ggplot(small)+geom_histogram(aes(x=price,fill=cut))
##position="dodge",将单一变量内部的比例分布分成不同柱画
ggplot(small)+geom_histogram(aes(x=price,fill=cut),position="dodge")

##position="fill"参数，均一化，突出变量内的比例分布
ggplot(small)+geom_histogram(aes(x=price,fill=cut),position="fill")

#######################################
##柱状图geom_bar,同样单变量，可以做出分布
##和直方图相比，柱状图适合对少量确定值进行统计
ggplot(small)+geom_bar(aes(x=clarity))

##通过stat参数可以调节柱状图的高度
ggplot()+geom_bar(aes(x=c(LETTERS[1:3]),y=1:3),stat="identity")
  
##############################            
##函数密度图 geom_density
##colour参数指定的是曲线的颜色
ggplot(small)+geom_density(aes(x=price,colour=cut))
##fill是往曲线下面填充颜色
ggplot(small)+geom_density(aes(x=price,fill=clarity))

###############################
## 箱线图 geom_boxplot
## 数据点比较少时，可以反应出比柱状图+errorbar更多信息

ggplot(small)+geom_boxplot(aes(x=cut,y=price,fill=color))


##画图就是在做映射，选择合适的几何对象（散点图、柱状图、直方图等）
##将不同的属性（例如xy轴数据，填充颜色等）写在不同几何对象括号中


##标尺,可以对属性中的值进行变换，例如取log10，或者换成彩虹色

ggplot(small)+geom_point(aes(x=carat,y=price,shape=cut,colour=color))+scale_y_log10()+scale_colour_manual(values=rainbow(7))

##统计变换，statistics,用法是stat_xxx, stat需要知道aes数据，所以将aes参数写在ggplot中
ggplot(small,aes(x=carat,y=price))+geom_point()+scale_y_log10()+stat_smooth()

##坐标系统，coordinante

##坐标翻转，coord_flip()
ggplot(small)+geom_bar(aes(x=cut,fill=cut))+coord_flip()

##转换为极坐标coord_polar()
ggplot(small)+geom_bar(aes(x=factor(1),fill=cut))+coord_polar(theta="y")

##靶心图
ggplot(small)+geom_bar(aes(x=factor(1),fill=cut))+coord_polar()


##风琴玫瑰图
ggplot(small)+geom_bar(aes(x=clarity,fill=cut))+coord_polar()


##图层，也是ggplot方便的原因.用+连接

##分面Facet，可以让我们对数据分组，一次对多组数据作图
ggplot(small,aes(x=carat,y=price))+geom_point(aes(colour=cut))+scale_y_log10()+facet_wrap(~cut)+stat_smooth()


##主题Theme
ggplot(small)+geom_boxplot(aes(x=cut,y=price,position="dodge",colour=color,fill=color))+theme_bw()+ggtitle("Price vs Cut")+xlab("Cut")+ylab("Price")
                           