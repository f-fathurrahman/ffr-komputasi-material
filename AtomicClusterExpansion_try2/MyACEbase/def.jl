
# Chris Rackaukas' macro generator
"""
Chris Rackaukas' macro generator

e.g.
```julia
@def mymac begin
   a = 2
end

@mymac
```

Used e.g. to define long lists of imports that are used in several
submodules.
"""
macro def(name, definition)
   return quote
      macro $(esc(name))()
         esc($(Expr(:quote, definition)))
      end
   end
end
