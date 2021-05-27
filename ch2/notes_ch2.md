## Some notes about chapter 2

### page 65:
What goes wrong if this step (i.e., $F1(x,\tau) = F1(x,\tau-1)$) 
is carried out in reverse order in $\tau$? 

If we would start at time $t$ with $\tau = 1$
and go to $\tau = \hat{\tau}$, the first assignment we would make would
be $F1(x,\tau=1) = F1(x,\tau=0)$, the next one $F1(x,\tau=2) = F1(x,\tau=1)$, 
where the latter value is $F1(x,\tau=0)$.

Hence, that would mean that the full array of F1 values gets assigned the
value $F1(x,\tau=0)$.

Conversely, when one goes backwards, the first assignment one makes is
$F1(x,\tau=\hat{\tau}) = F1(x,\tau=\hat{\tau}-1)$, the next one 
$F1(x,\tau=\hat{\tau} - 1) = F1(x,\tau=\hat{\tau} - 2)$. In this case,
you merely shift the contents of the vector 1 upwards.
