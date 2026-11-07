"""Separate Python implication search for minimal siphons (cross-check of the C++ enumeration)."""

from time import perf_counter

def enumerate_models(model, kind='siphon', progress=None):
    order=list(model.species);index={s:i for i,s in enumerate(order)}
    def mask(side):return sum(1<<index[s] for s,n in side.items() if n>0)
    pre=[mask(r.lhs) for r in model.reactions.values()];post=[mask(r.rhs) for r in model.reactions.values()]
    clauses=list(zip(post,pre)) if kind=='siphon' else list(zip(pre,post))
    result=[];states=0
    def closure(value):
        while True:
            old=value
            for t,r in clauses:
                if value&t and not value&r:value&=~t
            if value==old:return value
    for seed in range(len(order)):
        stack=[(1<<seed,(1<<seed)-1)]
        while stack:
            current,excluded=stack.pop();states+=1
            while True:
                before=(current,excluded)
                conflict=False
                for known in result:
                    if known&excluded:continue
                    missing=known&~current
                    if not missing:conflict=True;break
                    if missing&(missing-1)==0:excluded|=missing
                if conflict:break
                for trigger,required in clauses:
                    if required&current:continue
                    choices=required&~excluded
                    if not choices:excluded|=trigger
                    elif trigger&current and choices&(choices-1)==0:current|=choices
                if current&excluded:conflict=True;break
                if before==(current,excluded):break
            if conflict:continue
            choices=None
            for t,r in clauses:
                if t&current and not r&current:
                    available=r&~excluded
                    if choices is None or available.bit_count()<choices.bit_count():choices=available
            if choices is None:
                minimal=current;bits=current
                while bits:
                    p=bits&-bits;bits-=p
                    if minimal&p:
                        smaller=closure(minimal&~p)
                        if smaller:minimal=smaller
                result.append(minimal)
            else:
                while choices:
                    bit=choices&-choices;choices-=bit
                    stack.append((current|bit,excluded))
                    excluded|=bit
            if progress and states%10000==0:progress(states,len(result),len(stack))
    return [{order[i] for i in range(len(order)) if mask&(1<<i)} for mask in result],states

