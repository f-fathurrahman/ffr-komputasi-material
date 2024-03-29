
# --------------------------------------------------------------------------
# MyACE1.jl: Julia implementation of the Atomic Cluster Expansion
# Copyright (c) 2019 Christoph Ortner <christophortner0@gmail.com>
# Licensed under ASL - see ASL.md for terms and conditions.
# --------------------------------------------------------------------------



export PolyPairBasis


struct PolyPairBasis{TJ, NZ} <: IPBasis
   zlist::SZList{NZ}
   J::SMatrix{NZ, NZ, TJ}
   bidx0::SMatrix{NZ,NZ,Int}
end

fltype(pB::PolyPairBasis) = fltype(pB.J[1])   # assume they are all the same

Base.length(pB::PolyPairBasis{TJ, NZ}) where {TJ, NZ} = 
      sum( length(pB.J[i,j]) for i = 1:NZ for j = i:NZ )

# these make no sense - where are they used? 
Base.length(pB::PolyPairBasis, z::AtomicNumber, z0::AtomicNumber) = 
         length(pB, z2i(pB, z), z2i(pB, z0))

Base.length(pB::PolyPairBasis, iz::Integer, iz0::Integer) = 
         length(pB.J[iz0, iz])

_getJ(pB::PolyPairBasis, z::AtomicNumber, z0::AtomicNumber) = 
         _getJ(pB, z2i(pB, z), z2i(pB, z0))

_getJ(pB::PolyPairBasis, i::Integer, i0::Integer) = pB.J[i0, i]


zlist(pB::PolyPairBasis) = pB.zlist

function PolyPairBasis(J, species)
   @info "This PolyPairBasis constructor is called 40"
   return PolyPairBasis(J, ZList(species; static=true))
end

function PolyPairBasis(J::ScalarBasis, species::SZList)
   @info "This PolyPairBasis constructor is called 44"
   return PolyPairBasis( [J for i = 1:length(species), j = 1:length(species)], species)
end

function PolyPairBasis(J::AbstractMatrix, species::SZList)
   @info "This PolyPairBasis constructor is called 48"
   return PolyPairBasis(SMatrix{length(species), length(species)}(J), species)
end

function PolyPairBasis(J::SMatrix{NZ, NZ, TJ}, zlist::SZList{NZ}) where {NZ, TJ <: ScalarBasis}
   @info "This PolyPairBasis constructor is called 53"
   return PolyPairBasis(zlist, J, get_bidx0(J, zlist))
end

function get_bidx0(J::SMatrix, zlist::SZList{NZ}) where {NZ}
   bidx0 = fill(zero(Int), (NZ, NZ))
   i0 = 0
   for i = 1:NZ, j = i:NZ
      bidx0[i,j] = i0
      bidx0[j,i] = i0
      i0 += length(J[i, j])
      # just a quick sanity check while we are at it ... 
      @assert length(J[i,j]) == length(J[j,i])
   end
   return SMatrix{NZ, NZ, Int}(bidx0)
end

function scaling(pB::PolyPairBasis, p)
   ww = zeros(Float64, length(pB))
   for iz0 = 1:numz(pB), iz = iz0:numz(pB)
      idx0 = _Bidx0(pB, iz0, iz)
      for n = 1:length(pB.J[iz0, iz])
         # TODO: very crude, can we do better?
         #       -> need a proper H2-orthogonality?
         ww[idx0+n] = n^p
      end
   end
   return ww
end

==(B1::PolyPairBasis, B2::PolyPairBasis) =
      (B1.J == B2.J) && (B1.zlist == B2.zlist)

MyJuLIP.cutoff(pB::PolyPairBasis) = maximum(cutoff.(pB.J))

write_dict(pB::PolyPairBasis) = Dict(
      "__id__" => "MyACE1_PolyPairBasis",
          "Pr" => write_dict.(pB.J[:]),
       "zlist" => write_dict(pB.zlist) )

function read_dict(::Val{:MyACE1_PolyPairBasis}, D::Dict) 
   zlist = read_dict(D["zlist"])
   NZ = length(zlist)
   J = reshape(read_dict.(D["Pr"]), (NZ, NZ))
   return PolyPairBasis(J, zlist)
end

alloc_temp(pB::PolyPairBasis, args...) = (
              J = alloc_B.(pB.J),
          tmp_J = alloc_temp.(pB.J)  )

alloc_temp_d(pB::PolyPairBasis, args...) =  (
             J = alloc_B.( pB.J),
         tmp_J = alloc_temp.(pB.J),
            dJ = alloc_dB.(pB.J),
        tmpd_J = alloc_temp_d.(pB.J)  )

"""
compute the zeroth index of the basis corresponding to the potential between
two species zi, zj; as precomputed in `PolyPairBasis.bidx0`
"""
_Bidx0(pB, zi, zj) = pB.bidx0[ z2i(pB, zi), z2i(pB, zj) ]
_Bidx0(pB, i::Integer, j::Integer) = pB.bidx0[ i, j ]

function energy(pB::PolyPairBasis, at::Atoms{T}) where {T}
   E = zeros(T, length(pB))
   tmp = alloc_temp(pB)
   for (i, j, R) in pairs(at, cutoff(pB))
      r = norm(R)
      Zi, Zj = at.Z[i], at.Z[j]
      Ii = z2i(pB, Zi); Ij = z2i(pB, Zj)
      J = tmp.J[Ii, Ij]
      evaluate!(J, tmp.tmp_J[Ii, Ij], pB.J[Ii, Ij], r, Zi, Zj)
      idx0 = _Bidx0(pB, Zi, Zj)
      for n = 1:length(pB.J[Ii, Ij])
         E[idx0 + n] += 0.5 * J[n]
      end
   end
   return E
end

function forces(pB::PolyPairBasis, at::Atoms{T}) where {T}
   F = zeros(JVec{T}, length(at), length(pB))
   tmp = alloc_temp_d(pB)
   for (i, j, R) in pairs(at, cutoff(pB))
      r = norm(R)
      Zi, Zj = at.Z[i], at.Z[j]
      Ii = z2i(pB, Zi); Ij = z2i(pB, Zj)
      dJ = tmp.dJ[Ii, Ij]
      evaluate_d!(tmp.J[Ii, Ij], dJ, tmp.tmpd_J[Ii, Ij], pB.J[Ii, Ij], r, Zi, Zj)
      idx0 = _Bidx0(pB, Zi, Zj)
      for n = 1:length(pB.J[Ii, Ij])
         F[i, idx0 + n] += 0.5 * dJ[n] * (R/r)
         F[j, idx0 + n] -= 0.5 * dJ[n] * (R/r)
      end
   end
   return [ F[:, iB] for iB = 1:length(pB) ]
end

function virial(pB::PolyPairBasis, at::Atoms{T}) where {T}
   V = zeros(JMat{T}, length(pB))
   tmp = alloc_temp_d(pB)
   for (i, j, R) in pairs(at, cutoff(pB))
      r = norm(R)
      Zi, Zj = at.Z[i], at.Z[j]
      Ii = z2i(pB, Zi); Ij = z2i(pB, Zj)
      dJ = tmp.dJ[Ii, Ij]
      evaluate_d!(tmp.J[Ii, Ij], dJ, tmp.tmpd_J[Ii, Ij], pB.J[Ii, Ij], r, Zi, Zj)
      idx0 = _Bidx0(pB, Zi, Zj)
      for n = 1:length(pB.J[Ii, Ij])
         V[idx0 + n] -= 0.5 * (dJ[n]/r) * R * R'
      end
   end
   return V
end


# ----------- 
# an implementation needed to use pair potentials as site-based descriptors 

function evaluate!(Ei::Vector{Float64}, tmp, basis::PolyPairBasis, Rs, Zs, z0)
   for (rr, z) in zip(Rs, Zs)
      r = norm(rr) 
      Ii = z2i(basis, z0)
      Ij = z2i(basis, z)
      J = tmp.J[Ii, Ij]
      evaluate!(J, tmp.tmp_J[Ii, Ij], basis.J[Ii, Ij], r, z0, z)
      idx0 = _Bidx0(basis, z0, z)
      for n = 1:length(basis.J[Ii, Ij])
         Ei[idx0 + n] += 0.5 * J[n]
      end
   end
   return Ei 
end


function evaluate_d!(dB::Matrix{<: JVec}, tmp, basis::PolyPairBasis, Rs, Zs, z0)
   fill!(dB, zero(eltype(dB)))
   for (j, (rr, z)) in enumerate(zip(Rs, Zs))
      r = norm(rr) 
      Ii = z2i(basis, z0)
      Ij = z2i(basis, z)
      dJ = tmp.dJ[Ii, Ij]
      evaluate_d!(tmp.J[Ii, Ij], dJ, tmp.tmpd_J[Ii, Ij], basis.J[Ii, Ij], r, z0, z)
      idx0 = _Bidx0(basis, z0, z)
      for n = 1:length(dJ)
         dB[idx0 + n, j] += 0.5 * dJ[n] * (rr/r)  # += 0.5 * J[n]
      end
   end
   return dB  
end