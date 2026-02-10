#!/usr/bin/env Rscript
# SPDX-FileCopyrightText: 2025 Ovidio Garcia-Oliva
# SPDX-License-Identifier: CC-BY-4.0
# SPDX-FileContributor Ovidio Garcia-Oliva <ovidio.garcia@hereon.de
#
# This script creates contour plots for limnological variables 

library(rLakeAnalyzer)
library(lubridate)
library(plotly)
library(cmocean)

setwd(system('pwd',intern = T))
par(las=1)

citation('rLakeAnalyzer')
citation('plotly')



CWG.temp.all = read.csv("./data/Lake_Atitlan_water_quality_CWG_2010-2024.csv",comment.char = '#')
CWG.temp.all$date = ymd(CWG.temp.all$date)

###########

bg.color = 'white' 
g.color = '#EEEEEE' 
cont.col = 'black' 
the.type = "contour"
#the.type='heatmap'
month.list = c('January','February','March','April','May','June','July','Aμgust','September','October','November','December')
nutrient.color = list(c(0.1 * 0:9, 0.99,1), colorRampPalette(c("#FFD700", "#00B7C2", "#247BA0", "#3D348B"))(12))


scale.factor = 1
the.w = 1000*scale.factor
the.h = 200*scale.factor
the.scale = 2

t1 = list(
  family = "Helvetica",
  #size = 14,
  color = "black")

###########

max = function(...) base::max(...,na.rm = T)
min = function(...) base::min(...,na.rm = T)

add_vline = function(p, x, ...) {
  
  if(!is.null(p$x$layoutAttrs)){
    index = unname(which(sapply(p$x$layoutAttrs, function(x) 
      !is.null(x$shapes))))
  } else {
    index = integer()
  }
  
  l_shape = list(
    type = "line",
    y0 = 0, y1 = 1, yref = "paper", # i.e. y as a proportion of visible region
    x0 = x, x1 = x,
    line = list(
      ...
    ),
    layer = "above"
  )
  
  if(length(index) > 0){
    shapes = p$x$layoutAttrs[[index]]$shapes
    shapes[[length(shapes) + 1]] = l_shape
    p$x$layoutAttrs[[index]]$shapes = shapes
  } else {
    p = plotly::layout(
      p = p,
      shapes = list(l_shape)
    )
  }
  p
}

column.filler = function(db){
  col.z = NULL
  col.depth = NULL
  col.date = NULL
  for(i in unique(db$date)){
    column = db[db$date==i,]
    if(length(column$depth)<2)next
    #x = min(column$depth):max(column$depth)
    x = seq(round(min(column$depth)),round(max(column$depth)+0.1),1)
    y = approx(column$depth,column$z,x)$y
    col.depth = c(col.depth,x)
    col.z = c(col.z,y)
    col.date = c(col.date,rep(i,length(x)))
  }
  db.r = as.data.frame(list(date=as.Date(col.date,'1970-01-01'),depth=col.depth,z=col.z))
  return(db.r)
}

contour.ploter = function(dataset,contour.data,max.depth=300,color='YlOrRd',rev=T){
  #dataset = dataset[dataset$depth<=max.depth,]
  
  fig0 = plot_ly(z = dataset$z,
                 x = dataset$date,
                 y = dataset$depth, 
                 connectgaps = F,
                 type = the.type,
                 colorscale=color, reversescale=rev,
                 line = list(width = 0., color = "BLACK"),
                 autocontour = F,
                 line = list(smoothing = 1.3),
                 contours = contour.data
  ) 
  
  for(i in 2011:2024){
    fig0 = fig0 %>% 
      add_vline(paste0(i,'-01-01'),width=0.75,color='#F7F7F7')  
  }
  
  fig0 = fig0 %>% 
    add_vline(paste0(2010,'-01-01'),width=2,color='#444444')  
  
  fig0 = fig0 %>% 
    add_vline(paste0(2025,'-01-01'),width=2,color='#444444') 
  
  fig0 = fig0 %>% layout(font=t1,  margin = list(t = 30) )
  
  fig0 = fig0 %>%
    colorbar(
      bgcolor=alpha(cont.col,.0),
      bordercolor=cont.col,
      tickcolor=cont.col,
      #borderwidth=0,
      orientation='h',
      xanchor='left',
      yanchor='bottom',
      xref='paper',
      yref='paper',
      x=0,
      y=0,
      len=0.15*the.w,
      lenmode='pixels',
      thickness=0.05*the.h,
      thicknessmode='pixels',
      #outlinewidth=2,
      outlinecolor=cont.col,
      nticks=4,
      tickfont=list(family='Carlito',weight='bold',lineposition='over',color=cont.col)
    )
  
  fig0 = fig0 %>%
    layout(
      xaxis = list(title = "",range=as.Date(c('2010-01-01','2025-01-01')))
    )
  
  fig0 = fig0 %>%
    add_trace(
      x = as.Date(c('2010-01-01','2025-01-01')),
      y = 0*c(-1,-1)+300+0*max(dataset$depth),
      type = "scatter",
      mode = "lines",
      line = list(color = "#444444", width = 3,connectgaps = F,shape = "hvh"),
      name='botttom'
      #point = list(color='white')
      #yaxis = "y2"  # Use secondary y-axis
    )
  
  
  return(fig0)
}

################################################################################
## CWG temp profile

CWG.temp.all$z = CWG.temp.all$temperature
CWG.temp = column.filler(na.omit(subset(CWG.temp.all,select=c(date,depth,z))))

range(CWG.temp$z)
contour.data = list(
  start = 19,
  end = 27,
  size = 0.1,
  showlabels=F)

temp.color = list(0.1 * 0:10, rev(cmocean('thermal')(11)))
fig0 = contour.ploter(CWG.temp,contour.data,color=temp.color,max.depth = 300)

fig0 = fig0 %>%
  layout(
    title = list(text="a. water temperature (°C)",
                 x=0,
                 xanchor='left',
                 xref='paper',
                 y=1.001,
                 yanchor='bottom',
                 yref='paper'
    ),
    xaxis = list(title = "",gridcolor = g.color, showgrid=F),
    yaxis = list(title = "depth (m)", showgrid=F, gridcolor = g.color,range = c(300,0)),  # Reverse y-axis
    plot_bgcolor = bg.color  # Set background color
  )


fig0
save_image(fig0,'./out/CWG_temp.pdf',width = the.w, height = the.h,scale=the.scale)
htmlwidgets::saveWidget(widget=fig0,'./out/CWG_temp.html', selfcontained = T)

################################################################################
## CWG DO profile

CWG.temp.all$z = CWG.temp.all$dissolved_oxygen
CWG.temp = column.filler(na.omit(subset(CWG.temp.all,select=c(date,depth,z))))

contour.data = list(
  start = 0,
  end = 8.1,
  size = 0.2,
  showlabels=F)

#DO.col = list(c(0, 0.25,0.5,0.75, 1), c('gold', 'snow','lightgray','darkgray','darkred'))
#DO.col='RdBu'
DO.col='Jet'
#DO.col = list(c(0, 0.25,0.5,0.75, 1), rev(cmocean('oxy')(5)))
DO.col = list(c(0, 0.25, 0.5,0.75, 1), c('yellow','gold','lightgray','darkgray','darkred'))

fig0 = contour.ploter(CWG.temp,contour.data,max.depth = 300,color=DO.col,rev=T)

if(F){
  fig0 = fig0 %>%
    add_trace(
      x =CWG$date ,
      y = CWG$oxycline,
      type = "scatter",
      mode = "points",
      line = list(color = "black", width = 4,connectgaps = F,shape = "hvh"),
      name = 'oxicline'
      #point = list(color='white')
      #yaxis = "y2"  # Use secondary y-axis
    )
}

fig0 = fig0 %>%
  layout(
    title = list(text="b. dissolved oxygen (mg L<sup>-1</sup>)",
                 x=0,
                 xanchor='left',
                 xref='paper',
                 y=1.001,
                 yanchor='bottom',
                 yref='paper'
    ),
    xaxis = list(title = "",gridcolor = g.color, showgrid=F),
    yaxis = list(title = "depth (m)", showgrid=F, gridcolor = g.color,range = c(300,0)),  # Reverse y-axis
    plot_bgcolor = bg.color  # Set background color
  )

fig0
save_image(fig0,'./out/CWG_DO.pdf',width = the.w, height = the.h,scale=the.scale)
htmlwidgets::saveWidget(widget=fig0,'./out/CWG_DO.html', selfcontained = T)

#stop-here
################################################################################
## CWG PO4 profile

CWG.temp.all$z = CWG.temp.all$ortho.phosphate
CWG.temp = column.filler(na.omit(subset(CWG.temp.all,select=c(date,depth,z))))
CWG.temp$z = log10(CWG.temp$z)

contour.data = list(
  start = 0,
  end = max(CWG.temp$z),
  size = 0.1,
  showlabels=F)

fig0 = contour.ploter(CWG.temp,contour.data,max.depth = 300,color=nutrient.color,rev=T)

fig0 = fig0 %>%
  layout(
    title = list(text="e. soluble reactive phosphorus (SRP) [PO<sub>4</sub>] (μg-P L<sup>-1</sup>)",
                 x=0,
                 xanchor='left',
                 xref='paper',
                 y=1.001,
                 yanchor='bottom',
                 yref='paper'
    ),
    xaxis = list(title = "",gridcolor = g.color, showgrid=F),
    yaxis = list(title = "depth (m)", showgrid=F, gridcolor = g.color,range = c(300,0)),  # Reverse y-axis
    plot_bgcolor = bg.color  # Set background color
  )

fig0 = fig0%>%
  colorbar(
    tickprefix='10<sup>',
    ticksuffix='</sup>'
  )


fig0
save_image(fig0,'./out/CWG_PO4.pdf',width = the.w, height = the.h,scale=the.scale)
htmlwidgets::saveWidget(widget=fig0,'./out/CWG_PO4.html', selfcontained = T)
#stop-here
################################################################################
## CWG NH4 profile

CWG.temp.all$z = CWG.temp.all$ammonia
CWG.temp = column.filler(na.omit(subset(CWG.temp.all,select=c(date,depth,z))))
CWG.temp$z = log10(CWG.temp$z)

contour.data = list(
  start = 0,
  end = max(CWG.temp$z),
  size = 0.1,
  showlabels=F)

fig0 = contour.ploter(CWG.temp,contour.data,max.depth = 300,color=nutrient.color,rev=T)

fig0 = fig0 %>%
  layout(
    title = list(text="d. ammonia [NH<sub>4</sub>] (μg-N L<sup>-1</sup>)",
                 x=0,
                 xanchor='left',
                 xref='paper',
                 y=1.001,
                 yanchor='bottom',
                 yref='paper'
    ),
    xaxis = list(title = "",gridcolor = g.color, showgrid=F),
    yaxis = list(title = "depth (m)", showgrid=F, gridcolor = g.color,range = c(300,0)),  # Reverse y-axis
    plot_bgcolor = bg.color  # Set background color
  )

fig0 = fig0%>%
  colorbar(
    tickprefix='10<sup>',
    ticksuffix='</sup>'
  )


fig0
save_image(fig0,'./out/CWG_NH4.pdf',width = the.w, height = the.h,scale=the.scale)
htmlwidgets::saveWidget(widget=fig0,'./out/CWG_NH4.html', selfcontained = T)

################################################################################
## CWG NOX profile

CWG.temp.all$z = CWG.temp.all$nitrate_plus_nitrite
CWG.temp = column.filler(na.omit(subset(CWG.temp.all,select=c(date,depth,z))))
CWG.temp$z = log10(CWG.temp$z)

contour.data = list(
  start = 0,
  end = max(CWG.temp$z),
  size = 0.1,
  showlabels=F)

fig0 = contour.ploter(CWG.temp,contour.data,max.depth = 300,color=nutrient.color,rev=T)

fig0 = fig0 %>%
  layout(
    title = list(text="c. nitrate + nitrite [NO<sub>3</sub>+NO<sub>2</sub>] (μg-N L<sup>-1</sup>)",
                 x=0,
                 xanchor='left',
                 xref='paper',
                 y=1.001,
                 yanchor='bottom',
                 yref='paper'
    ),
    xaxis = list(title = "",gridcolor = g.color, showgrid=F),
    yaxis = list(title = "depth (m)", showgrid=F, gridcolor = g.color,range = c(300,0)),  # Reverse y-axis
    plot_bgcolor = bg.color  # Set background color
  )

fig0 = fig0%>%
  colorbar(
    tickprefix='10<sup>',
    ticksuffix='</sup>'
  )


fig0
save_image(fig0,'./out/CWG_NOx.pdf',width = the.w, height = the.h,scale=the.scale)
htmlwidgets::saveWidget(widget=fig0,'./out/CWG_NOx.html', selfcontained = T)

################################################################################
## CWG TP profile

CWG.temp.all$z = CWG.temp.all$total_phosphorus
CWG.temp = column.filler(na.omit(subset(CWG.temp.all,select=c(date,depth,z))))
CWG.temp$z = log10(CWG.temp$z)

contour.data = list(
  start = 0,
  end = max(CWG.temp$z),
  size = 0.1,
  showlabels=F)

fig0 = contour.ploter(CWG.temp,contour.data,max.depth = 300,color=nutrient.color,rev=T)

fig0 = fig0 %>%
  layout(
    title = list(text= "f. total phosphorus [TP] (μg-P L<sup>-1</sup>)",
                 x=0,
                 xanchor='left',
                 xref='paper',
                 y=1.001,
                 yanchor='bottom',
                 yref='paper'
    ),
    xaxis = list(title = "",gridcolor = g.color, showgrid=F),
    yaxis = list(title = "depth (m)", showgrid=F, gridcolor = g.color,range = c(300,0)),  # Reverse y-axis
    plot_bgcolor = bg.color  # Set background color
  )

fig0 = fig0%>%
  colorbar(
    tickprefix='10<sup>',
    ticksuffix='</sup>'
  )


fig0
save_image(fig0,'./out/CWG_TP.pdf',width = the.w, height = the.h,scale=the.scale)
htmlwidgets::saveWidget(widget=fig0,'./out/CWG_TP.html', selfcontained = T)
