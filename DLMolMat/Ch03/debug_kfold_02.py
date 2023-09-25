N = 110
N = 110
for k in range(2, 25):
    print("Number of folds: ", k)
    splits = list(range(0, N + N // k, N // k))
    print("splits = ", splits)
    #for i in range(k):
    #    print("start splits[i] = ", splits[i], " splits[i+1] = ", splits[i+1])

