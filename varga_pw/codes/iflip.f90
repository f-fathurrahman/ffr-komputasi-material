function iflip(i,n)
  IMPLICIT NONE 
  INTEGER  :: i,n,iflip

  if(i.gt.n/2) then
    iflip=i-n-1
  else
    iflip=i-1
  endif
end function iflip

function iflip_inv(i,n)
  IMPLICIT NONE 
  INTEGER  :: i,n,iflip_inv

  iflip_inv=i+1
  if(i.lt.0) iflip_inv=iflip_inv+n

end function iflip_inv



