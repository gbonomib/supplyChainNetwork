library(igraph)
library(dplyr)
library(MCMCpack)

nTiers = 5
nTiers = nTiers + 2

tiers = vector('list', nTiers)


for(i in (seq_len(nTiers))){
  
  tiers[[i]] = seq_len(rbinom(1,10,0.6) + 1) 
  
}

tiers[1] = 1
tiers[length(tiers)] = 1

for(i in seq_len(nTiers)[-1]){
  
  tiers[[i]] = tiers[[i]] + max(tiers[[i-1]])
  
}


i = 1
j = 1

# connectionThreshold = 0.3
edgelist = c()

for(i in seq_len(nTiers)[-nTiers]){
  
  # connectionThreshold = ifelse(i != 1 & i != length(tiers), 0.3, 1)
  
  for(j in tiers[[i]]){
    
    # links = vector('logical', length(tiers[[i+1]]))
    # 
    # while(sum(links) == F){
    #   
    #   for(k in sample(seq_along(tiers[[i+1]]))){
    #     
    #     links[k] = ifelse(runif(1,0,1) > 1 - connectionThreshold, k + max(tiers[[i]]), FALSE)  
    #     
    #   }
    # }
    # 
    # links = links[links!=0]
    
    weights = vector('logical',length(tiers[[i+1]]))
    
    if(i == 1){
      
      while(prod(weights) == 0){
        
        weights = trunc(rdirichlet(1,runif(length(weights),0,1))*1e2)/1e2
          
        # for(k in sample(seq_along(weights))){              
        #   weights[k] = min(runif(1,0,1),
        #                    max(0, 1 - sum(weights, na.rm = T)))
        #   
        # }
      }
    }else{
      
      weights = trunc(rdirichlet(1,runif(length(weights),0,1))*1e2)/1e2
      
      # for(k in sample(seq_along(weights))){
      #   
      #   weights[k] = min(runif(1,0,1),
      #                    max(0, 1 - sum(weights, na.rm = T)))
      #  
      #}
      
    }
    
    weights = weights %>% c()
    weights = weights/sum(weights)*1e2/1e2
    
    edgelistTmp = cbind(rep(j,length(tiers[[i+1]])),tiers[[i+1]],weights)
    
    edgelist = rbind(edgelist,edgelistTmp)
    
    
    
  }}


supplyChain = graph_from_data_frame(edgelist)
plot(supplyChain)

E(supplyChain)$weights
V(supplyChain)$name

# at this point the supply chain is istantiated, the connection between the nodes created and weighted
# so that each company allocates 100% of its production downstream

# the size and the production capacity of each company should be proportional to the
# weights that such company receives in input (something like sum(upstream_edgesWeights))

V(supplyChain)$strength = strength(supplyChain, vids =V(supplyChain),
                                   mode = 'in',
                                   weights = E(supplyChain)$weights)

# normalize weights

V(supplyChain)$strength[1] = 1 # this needs to be set equal to one
V(supplyChain)$size[1] = 1 # this needs to be set equal to one

for(company in V(supplyChain)[-1]){
  
  suppliers = neighbors(supplyChain, company, mode = 'in')
  
  suppliersStrength = get.vertex.attribute(supplyChain, 'size', suppliers)
  
  connections = incident_edges(supplyChain, company, mode = "in") %>% unlist()
  
  connectionsStrength = get.edge.attribute(supplyChain, "weights", connections)
  
  V(supplyChain)$size[company] = sum(suppliersStrength * connectionsStrength)
  
}

plot(supplyChain)

for(company in V(supplyChain)[-last(V(supplyChain))]){
  
  connections = incident_edges(supplyChain, company, mode = 'out') %>% unlist()
  
  E(supplyChain)$flow[connections] = E(supplyChain)$weights[connections] * V(supplyChain)$size[company]
  
}  

#TODO check that the size is correct. i presume that the size tierwise shuold be always the same but i am not sure that it is

V(supplyChain)$size[1] = 1e-5
V(supplyChain)$size[last(V(supplyChain))] = 1e-5


supplyChain = delete_vertices(supplyChain, V(supplyChain)[V(supplyChain)$size == 0]) 

plot(supplyChain,
     vertex.size = V(supplyChain)$size*30,
     edge.arrow.size = 0.1,
     edge.width = E(supplyChain)$flow*30)

### now we have a overall plot of the supply chain, inclusding the relative size of the companies and
### the edges are weighted according to the flow which pass by them

### set max and min demand forecast

# dailyPieces = 1000

annualDemand = 365*1e3
demandSd = 100
maxDemand = 365*1e3
minDemand = 365*1e3*0.8

annualDemand = rnorm(1, max(min(annualDemand,maxDemand), minDemand), 100) # bounded markowchain setup

dailyPieces = annualDemand/365

### remember you will need to istantiate the poisson process for the realised demand
### the order is a R object (list?) which contains size and id

orderIds = 0

ordersNumber = dailyPieces
orderIds = last(orderIds) + seq_len(ordersNumber)

orders = data.frame('orderId' = orderIds, 'orderSize' = runif(ordersNumber,0,1))


### the attributes of each and every node should fall in one of two categories: OPERATIONAL or FINANCIAL

### set max productive capacity accoring to size.
### this should be equal to the max daily demand times the size of the company

V(supplyChain)$maxProdCapacity = ceiling(maxDemand/365*V(supplyChain)$size)
V(supplyChain)$prodCapacity = V(supplyChain)$maxProdCapacity
V(supplyChain)$stress = V(supplyChain)$maxProdCapacity*0.01

### set initial wealth according to size

V(supplyChain)$wealth = 1e6*V(supplyChain)$size

### set receivable days for every customer (according to the size of the specitic customer?)
### this may be initially set equal for everyone

V(supplyChain)$DSO = 40

### istantiate the attribute 'queue' for kepping non produced items at stock
### set max queue level according to size

V(supplyChain)$maxInboundQueue = ceiling(V(supplyChain)$size*1000)
V(supplyChain)$inboundQueueFree = V(supplyChain)$maxInboundQueue
V(supplyChain)$inboundQueue = list(data.frame('orderId' = integer(), 'orderSize' = numeric()))

V(supplyChain)$maxOutboundQueue = ceiling(V(supplyChain)$size*1000)
V(supplyChain)$outboundQueueFree = V(supplyChain)$maxOutboundQueue
V(supplyChain)$outboundQueue = list(data.frame('orderId' = integer(), 'orderSize' = numeric()))

### istantiate the attribute 'work in progess' for keeping the wip (in this stage the production is assumed sequential with a FIFO logic therefore only one item will fall in this bucket)

V(supplyChain)$maxWip = 1
V(supplyChain)$wip = list(data.frame('orderId' = integer(), 'orderSize' = numeric()))

### instantiate list of orders to be passed upstream

V(supplyChain)$ordersId = list(data.frame('orderId' = integer()))
V(supplyChain)$orders = list(data.frame('orderId' = integer(), 'orderSize' = numeric()))
V(supplyChain)[V(supplyChain)$name == length(V(supplyChain))]$orders[[1]] = orders
V(supplyChain)[V(supplyChain)$name == length(V(supplyChain))]$ordersId[[1]] = orderIds

V(supplyChain)$ordersIssued = list(data.frame('orderId' = integer(), 'orderSize' = numeric()))

### istantiate list of payables. order ids will be kept here untili they are paid in full
V(supplyChain)$payables = list(data.frame('orderId' = integer(), 'orderSize' = numeric(), "amountDue" = integer(), "DSO" = integer()))
V(supplyChain)$price = 1
V(supplyChain)$DPO = 30

### istantiate list of receivables. order ids will be kept here until they are paid in full

### istantiate a variable for storing the oders to be produced in the following period and to whom should they go (order id + customer name)
### this is updated recursively from downstream MISSA DI NO!

### initialize the queue, or set the maximum queue level
V(supplyChain)[V(supplyChain)$name == 37]$inboundQueue[[1]] = orders #temporary initialization, just for testing


### set stress level (i.e. variability of productive capacity) for everyone
### set specific distress level

### set cash level (cash principle)
V(supplyChain)$netCash = 1000

### set profit level (accrual principle)
V(supplyChain)$pnl = 1000

### bucket time in discrete periods

periods = 2500

V(supplyChain)$size[1] = 1
V(supplyChain)$size[last(V(supplyChain))] = 1

### THIS IS IMPORTANT: here must go the initial iteration where only the ORDER LOOP is run to populate the orders
### keep in mind that in the first n loops (n = number of tiers in the supplyChain) not every company will produce
### something, because no pieces have arrived yet to produce. SO PERHAPS THE ORDER LOOP MIGHT BE THE FIRST LOOP INSTEAD OF BEING THE 3rd! check that possibility


## eventually the all the dynamics of the model should evolve inside this for loop

for(i in periods){
  
  # HERE update demand from the market and orders from the markets with new id. The markets always buy what it has ordered (this assumption may be even released later)
  
  ## what happens to each and every company in a single period shall be written here
  ## there shall be multiple for loops: in the first one everythin is produced, in the second one purchases are issued in the third one sales are done in the fourth one payments happen
  ## transaction of goods and money and everythin else happen. this is better to ensure generalizability
  ## here the assumpion is that what arrives in input at the period t will not be produced that day.
  ## if this assumption can be released it is better but i do not se how it would be possible now.
  
  # 3. ORDER LOOP.
  # TODO clearly distinguish between to be ordered and already ordered.
  # this distinction is necessary to manage the ORDER and the SALE loop (i.e. where do i retrieve the ordered pieces from upstream?)
  
  for(company in sort(V(supplyChain), decreasing = T)[-length(V(supplyChain))]){ # everything should be pulled from downstream, so the loop should start from the last node
    
    upstream = neighbors(supplyChain, V(supplyChain)[V(supplyChain)$name == company], mode = 'in')
    
    # if(upstream == 1){break}
    
    ordersId = V(supplyChain)[V(supplyChain)$name == company]$ordersId[[1]] = V(supplyChain)[V(supplyChain)$name == company]$orders[[1]]$orderId
    
    cumulatedFlow = 1
    
    for(supplier in upstream){
      
      orders = V(supplyChain)[V(supplyChain)$name == company]$orders[[1]]
      
      flow = E(supplyChain)[V(supplyChain)[V(supplyChain)$name == company] %--% V(supplyChain)[V(supplyChain)$name == supplier]]$flow
      
      flow = flow/V(supplyChain)[V(supplyChain)$name == company]$size
      
      x = x + flow
      
      ordered = sample_frac(orders,min(1,flow/cumulatedFlow))
      
      V(supplyChain)[V(supplyChain)$name == supplier]$orders[[1]] = 
        rbind(V(supplyChain)[V(supplyChain)$name == supplier]$orders[[1]],
              ordered)
      
      cumulatedFlow = cumulatedFlow - flow
      
      V(supplyChain)[V(supplyChain)$name == company]$orders[[1]] = 
        V(supplyChain)[V(supplyChain)$name == company]$orders[[1]] %>% filter(!orderId %in% ordered$orderId)
      
      V(supplyChain)[V(supplyChain)$name == company]$ordersIssued[[1]] =
        rbind(V(supplyChain)[V(supplyChain)$name == company]$ordersIssued[[1]],
              ordered)
      
      
      # [-as.numeric(rownames(V(supplyChain)[V(supplyChain)$name == supplier]$orders[[1]])),]
      
    }
  } # this comes to an end without breaking
  
  # 1. PRODUCTION LOOP
    # add a brief description
  
  for(company in V(supplyChain)[-c(1,length(V(supplyChain)))]){ # everyone produces something except for the first and the last node
    
    # update of production capacity
    V(supplyChain)[V(supplyChain)$name == company]$prodCapacity = 
      min(rnorm(1,
                V(supplyChain)[V(supplyChain)$name == company]$prodCapacity,
                V(supplyChain)[V(supplyChain)$name == company]$stress),
          V(supplyChain)[V(supplyChain)$name == company]$maxProdCapacity)
    
    # production of what is in the queue
    availableProductionCapacity =
      V(supplyChain)[V(supplyChain)$name == company]$prodCapacity
    
    inboundQueue = V(supplyChain)[V(supplyChain)$name == company]$inboundQueue[[1]]$orderId
    
    if(prod(!is.na(inboundQueue)) == 1){
      
      for(wip in inboundQueue){
        
        inProcess = V(supplyChain)[V(supplyChain)$name == company]$inboundQueue[[1]] %>% filter(orderId == wip)
        
        if(inProcess$orderSize - availableProductionCapacity < 0){
          
          availableProductionCapacity = availableProductionCapacity - inProcess$orderSize
          
          if(nrow(V(supplyChain)[V(supplyChain)$name == company]$outboundQueue[[1]]) <= V(supplyChain)[V(supplyChain)$name == company]$maxOutboundQueue){
            
            V(supplyChain)[V(supplyChain)$name == company]$outboundQueue[[1]] = 
              V(supplyChain)[V(supplyChain)$name == company]$outboundQueue[[1]] %>% rbind(inProcess)  # put produced pieces in the outbound list. there which order goes to whom?
            
            V(supplyChain)[V(supplyChain)$name == company]$inboundQueue[[1]] = 
              V(supplyChain)[V(supplyChain)$name == company]$inboundQueue[[1]] %>% filter(orderId!=inProcess$orderId)
          }else{
            
            break
            
          }
          
        }else{
          # V(supplyChain)[V(supplyChain)$name == company]$wip = wip - availableProductionCapacity
          break # let's for now assume that if a process cannot be fully produced in the period, it is not produced in that period
        } # TODO maybe this will slighly change considering that the object "order" has two values: id and size
        
      }
      
    }
  } ## this comes to an end without breaking
  
  # 2. SALE LOOP
    #TODO remember that in the first iteration an ORDER LOOP must be istantiated OR orders should be pre inserted in each company
    # empty each company outboundQueue and migrate the produced pieces in the inboundQueue downstream taking into account how much space there is
    # update payables in the downstream companies
  
  for(company in V(supplyChain)){
    
    downstream = neighbors(supplyChain, V(supplyChain)[V(supplyChain)$name == company], mode = 'out')
    
    # somehow the order i cannot get because not yet produced must be kept in the order list
    
    for(buyer in downstream){
      
      oldQueue = V(supplyChain)[V(supplyChain)$name == buyer]$inboundQueue[[1]]
      queueFree = V(supplyChain)[V(supplyChain)$name == buyer]$inboundQueueFree = 
        V(supplyChain)[V(supplyChain)$name == buyer]$maxInboundQueue - nrow(oldQueue)
      
      orders = V(supplyChain)[V(supplyChain)$name == buyer]$ordersIssued[[1]]
      
      newOrders = V(supplyChain)[V(supplyChain)$name == company]$outboundQueue[[1]] %>% filter(orderId %in% orders)
      
      # newOrders = V(supplyChain)[V(supplyChain)$name == company]$queue[[1]][
      #   V(supplyChain)[V(supplyChain)$name == company]$queue[[1]]$orderId =
      #     V(supplyChain)[V(supplyChain)$name == buyer]$orderId]
      
      if(nrow(newOrders) >= V(supplyChain)[V(supplyChain)$name == buyer]$inboundQueueFree){
        
        V(supplyChain)[V(supplyChain)$name == buyer]$inboundQueue[[1]] = rbind(oldQueue,
                                                                        newOrders)

        
      }else{
        
        # newOrders = sample_n(newOrders,queueFree) #this is a random sampling but ideally orders should be sampled in order
        
        newOrders =  head(newOrders,queueFree) # take the first orders to be processed
        
        V(supplyChain)[V(supplyChain)$name == buyer]$inboundQueue[[1]] = rbind(oldQueue,
                                                                               newOrders)
        
        # V(supplyChain)[V(supplyChain)$name == company]$outboundQueue[[1]] =
        #   V(supplyChain)[V(supplyChain)$name == company]$outboundQueue[[1]] %>% filter(orderId != newOrders$orderId)
 
      }
      
      V(supplyChain)[V(supplyChain)$name == company]$outboundQueue[[1]] =
        V(supplyChain)[V(supplyChain)$name == company]$outboundQueue[[1]] %>% filter(orderId != newOrders$orderId)
      
      
      V(supplyChain)[V(supplyChain)$name == buyer]$ordersIssued[[1]] =
        V(supplyChain)[V(supplyChain)$name == buyer]$ordersIssued[[1]] %>% filter(orderId != newOrders$orderId)
      
      ## somehow update the payable list!
      
      price = V(supplyChain)[V(supplyChain)$name == buyer]$price # do we let the buyer to set the price? ok for now...
      
      newOrders = newOrders %>% mutate(amountDue = orderSize*price, daysFromPurchase = 1)
      
      oldPayables = V(supplyChain)[V(supplyChain)$name == buyer]$payables[[1]]
      
      V(supplyChain)[V(supplyChain)$name == buyer]$payables[[1]] = rbind(oldPayables,
                                                                         newOrders)
      
    }
    
    # this instruction i do not understand
    # V(supplyChain)[V(supplyChain)$name == company]$ordersId[[1]] = V(supplyChain)[V(supplyChain)$name == company]$orders[[1]]$orderId
    

    
  } # this comes to an end without breaking
  

  
  # 4. PAYMENTS LOOP. (B2C is to be handled separately before the loop)
  
  ordersNumber = rnorm(1, dailyPieces)
  orderIds = last(orderIds) + seq_len(ordersNumber)
  
  orders = data.frame('orderId' = orderIds, 'orderSize' = runif(ordersNumber,0,1))
  marketPurchases = rpois(1, dailyPieces)
  
### BEGIN market purchase

  market = last(V(supplyChain))$name
  
  upstream = neighbors(supplyChain, V(supplyChain)[V(supplyChain)$name == market], mode = 'in')
  
  # V(supplyChain)[V(supplyChain)$name == company]$ordersId[[1]] =
  #   V(supplyChain)[V(supplyChain)$name == company]$ordersIssued[[1]]$orderId
  
  for(supplier in upstream){
    
    flow = E(supplyChain)[V(supplyChain)[V(supplyChain)$name == supplier] %--% V(supplyChain)[V(supplyChain)$name == market]]$flow

    nPurchases = ceiling(marketPurchases*flow)
    
    payables = head(V(supplyChain)[V(supplyChain)$name == supplier]$payables, nPurchases)
    
      for(id in payables$orderId){ ## istantiate payment loop
        
        toBePaid = payables$orderId %>% filter(orderId == id)
        
          ### here update cash levels
          
          V(supplyChain)[V(supplyChain)$name == supplier]$netCash[[1]] =
            V(supplyChain)[V(supplyChain)$name == supplier]$netCash[[1]] + toBePaid$amountDue
          
          V(supplyChain)[V(supplyChain)$name == company]$payables[[1]] = 
            V(supplyChain)[V(supplyChain)$name == company]$payables[[1]] %>% filter(!orderId == payables$orderId)
        
      }
    
    ### update DSO
    
    V(supplyChain)[V(supplyChain)$name == company]$payables[[1]]$DSO = 
      V(supplyChain)[V(supplyChain)$name == company]$payables[[1]]$DSO + 1
    
  }  
  
### END market purchases
  
  for(company in sort(V(supplyChain), decreasing = T)[-1]){ ## the "market" not is handled separately before this loop
    
    ### all payments after DPO should be settled here. To keep track of the days passed from the purchase a counter may be updated at each period
    ### the market pays when receives the product
    ### the market may not always purchase the good. indeed, orders depend from the annual forecast while purchases depend on the poisson process
    
    upstream = neighbors(supplyChain, V(supplyChain)[V(supplyChain)$name == company], mode = 'in')
    
    # V(supplyChain)[V(supplyChain)$name == company]$ordersId[[1]] =
    #   V(supplyChain)[V(supplyChain)$name == company]$ordersIssued[[1]]$orderId
    
    for(supplier in upstream){
      
      ### which orders must be paid to the given supplier. look for corresponfding ids in the payable list
      
      payables = V(supplyChain)[V(supplyChain)$name == company]$payables[[1]] %>%
        filter(DSO > V(supplyChain)[V(supplyChain)$name == company]$DPO)
      
      if(nrow(payables) > 0){
        
        for(id in payables$orderId){ ## istantiate payment loop
          
          ### here settle payables. keep in mind that order shall be paid from the lowest id only if there is enough cash
           ### also, even though one might choose to pay less expensive orders in case of cash shortage, this is not the case for this simulation, not at this point at least
        
          toBePaid = payables$orderId %>% filter(orderId == id)
          
          if(toBePaid$amountDue < V(supplyChain)[V(supplyChain)$name == supplier]$netCash[[1]]){
            
            ### here update cash levels
            
            V(supplyChain)[V(supplyChain)$name == supplier]$netCash[[1]] =
              V(supplyChain)[V(supplyChain)$name == supplier]$netCash[[1]] + toBePaid$amountDue
            
            V(supplyChain)[V(supplyChain)$name == company]$netCash[[1]] =
              V(supplyChain)[V(supplyChain)$name == company]$netCash[[1]] - toBePaid$amountDue
            
            V(supplyChain)[V(supplyChain)$name == company]$payables[[1]] = 
              V(supplyChain)[V(supplyChain)$name == company]$payables[[1]] %>% filter(!orderId == payables$orderId)
            
          }else{
            break
          }
        }
      }
      
      ### update DSO
      
      V(supplyChain)[V(supplyChain)$name == company]$payables[[1]]$DSO = 
        V(supplyChain)[V(supplyChain)$name == company]$payables[[1]]$DSO + 1
      
    }
    
    
  }
  
  
  
  # (which is the market) the market behaves as a poisson process whose average at t is the 
  # market demand forecast
  
  # update of the queue (are some orders to be refused? should a penalty be paid in shuch situation?)
  # this implies checking the outbound list of every supplier to see what is ready to be shipped.
  # do the company have enough cash to purchase what they have ordered?
  # contextually the the payable list should be updated.
  
  # update profit (purchases and sales) according to what as been sold/purchased
  # update cash (inflows and outflofs) according to cash movements
  # update DPD situation in the payable list. if an order defaults then it could perhaps be eliminated from this list
  
  # pass over to the supplier the orders requested for the following period (or more in general for the n-th period ahead? maybe not now)
  # according to the order received from downtream (whether it be the market itself or the business partners)
  
  #}
  
  # at the end of each and every period a summary table should be updated with the statistics of all
  # of the companies (or at least with those of the one of interest). the objective is to use those
  # later in the building of the predictive model and for plotting reasons
  
} 


####

xx = c()
yy = c()

for(i in V(supplyChain)){
  
  xx[i] = V(supplyChain)[V(supplyChain)$name == i]$ordersIssued[[1]] %>% nrow()
  yy[i] =  V(supplyChain)[V(supplyChain)$name == i]$orders[[1]] %>% nrow()
}

sum(xx[34:41])
sum(xx[25:33])
sum(xx[17:24])
sum(xx[12:16])

plot(xx)
plot(yy)
V(supplyChain)[V(supplyChain)$name == 11]$orders
