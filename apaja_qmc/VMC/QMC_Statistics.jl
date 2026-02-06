#
# Routines for collecting measurement samples from QMC simultion
# 
# Initialization: For a single measured value,
#    stat = init_stat(1, blocksize)
# and for a measured array of ndata values,
#    stat = init_stat(ndata, blocksize) 
# When stat.finished is true, and one block is done, you can
# pick up the QMC error estimate with
#     average, std, input_σ2, Nblocks  = get_stats(stat) 
# where input_σ2 is the variance^2 of the raw input data
#

module QMC_Statistics

export init_stat, add_sample!, get_stats, t_Stat

mutable struct t_StatData
    n    ::Int64  
    data ::Vector{Float64}
    data2 ::Vector{Float64}
    input_σ2 :: Float64
end

mutable struct t_Stat    
    nblocks   ::Int64
    blocksize ::Int64
    finished  ::Bool
    sample    ::t_StatData
    datablock :: Vector{t_StatData}
    t_Stat() = new()
end



function init_stat(datasize ::Int64, blocksize ::Int64; numblocks::Int64=100)
    stat = t_Stat()
    stat.blocksize = blocksize
    stat.nblocks = 0    
    stat.sample = t_StatData(0, zeros(datasize), zeros(datasize), 0.0)
    # data blocks    
    stat.datablock = Vector{t_StatData}(undef, numblocks) # was []
    stat.finished = false
    return stat
end

function add_sample!(stat ::t_Stat, dat)
    stat.finished  = false

    @. stat.sample.data  += dat
    @. stat.sample.data2 += dat^2
    stat.sample.n += 1
    
    if stat.sample.n == stat.blocksize
        # one full block collected, move average to block data
        # input data σ^2:
        N = stat.sample.n
        d2 = sum(stat.sample.data2)/N
        d  = sum(stat.sample.data)/N
        input_σ2 = d2-d^2
        #
        stat.nblocks += 1
        if stat.nblocks > length(stat.datablock)
            old_size =  length(stat.datablock)
            resize!(stat.datablock, old_size + 100) # unintialized elements in the end
        end 
        stat.datablock[stat.nblocks] = t_StatData(0, stat.sample.data./N, zeros(length(stat.sample.data)), input_σ2)
        # old, see "was" in init_stat
        #push!(stat.datablock,t_StatData(0, stat.sample.data./N, zeros(length(stat.sample.data)), input_σ2))        
        @. stat.sample.data = 0
        @. stat.sample.data2 = 0
        stat.sample.n = 0
        stat.finished = true
    end  
end



function get_stats(stat ::t_Stat)
    N = stat.nblocks
    if N==0
        println("get_stats: no data")
        return 0, 0, 0, 0
    end
    if length(stat.datablock[1].data)==1
        ave_1, std_1, input_σ2_1, N_1 = get_stats_1(stat ::t_Stat)
        return ave_1, std_1, input_σ2_1, N_1,  stat.datablock[N].data[1] 
    end
    ave = similar(stat.datablock[1].data) 
    ave2 = similar(ave)
    ave .= 0
    ave2 .= 0
    input_σ2 = 0.0
    @simd for i = 1:N
        d = stat.datablock[i].data
        @. ave += d
        @. ave2 += d^2
        input_σ2 += stat.datablock[i].input_σ2 
    end    
    @. ave /= N
    @. ave2 /= N
    input_σ2 /= N
    var = copy(ave)
    var2 = copy(ave)
    std = copy(ave)
    @. var2 = abs(ave2 - ave^2)
    @. var = sqrt(var2)
    @. std = var/sqrt(N)    
    return ave, std, input_σ2, N, stat.datablock[N].data 
end

function get_stats_1(stat ::t_Stat)
    N = stat.nblocks
    if N==0
        println("get_stats: no data")
        return nothing
    end
    ave = 0.0 
    ave2 = 0.0
    input_σ2 = 0.0 
    for i = 1:N
        d = stat.datablock[i].data[1]
        ave += d
        ave2 += d^2
        input_σ2 += stat.datablock[i].input_σ2 
    end 
    ave /= N
    ave2 /= N
    input_σ2 /= N
    var2 = abs(ave2 - ave^2)
    var = sqrt(var2) 
    std = var/sqrt(N)
    return ave, std, input_σ2, N
end


# not used
function extract_block_averages(stat::t_Stat)
    """get raw data as averages over each block, for all blocks so far"""
    N = stat.nblocks
    M = length(stat.datablock[1].data)
    blocks = Matrix{Float64}(undef, M, N)  # Each column = one block average
    for i in 1:N
        blocks[:, i] .= stat.datablock[i].data
    end
    return blocks
end


end
