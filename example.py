from mc2err import ExpValue

bernoulli_dist = ExpValue()

fake_data = [0,1]*50

bernoulli_dist.add(fake_data)

mean = bernoulli_dist.mean()
error = bernoulli_dist.error()
size = bernoulli_dist.size()
print("1st estimate:",mean,"+/-",error,"with",size,"effective data points")

more_fake_data = [0,1]*100

bernoulli_dist.add(more_fake_data)

mean = bernoulli_dist.mean()
error = bernoulli_dist.error()
size = bernoulli_dist.size()
print("2nd estimate:",mean,"+/-",error,"with",size,"effective data points")
