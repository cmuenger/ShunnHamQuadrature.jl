using StaticArrays

#Quadrature rulesfrom Shunn&Ham Paper

const pentatope1a = [
    SVector{5}([0.2 0.2 0.2 0.2 0.2] ),
]

const pentatope1w = [
    1.0
]



const pentatope5a = [
    SVector{5}([0.118350341907227 0.118350341907227 0.118350341907227 0.118350341907227 0.526598632371091] ),
    SVector{5}([0.118350341907227 0.118350341907227 0.118350341907227 0.526598632371091 0.118350341907227] ),
    SVector{5}([0.118350341907227 0.118350341907227 0.526598632371091 0.118350341907227 0.118350341907227] ),
    SVector{5}([0.118350341907227 0.526598632371091 0.118350341907227 0.118350341907227 0.118350341907227] ),
    SVector{5}([0.526598632371091 0.118350341907227 0.118350341907227 0.118350341907227 0.118350341907227] )
]

const pentatope5w = [
    0.2
    0.2
    0.2
    0.2
    0.2
]


const pentatope15a = [
    SVector{5}([0.0566663810400515 0.0566663810400515 0.0566663810400515 0.0566663810400515 0.773334475839794]), 
    SVector{5}([0.0566663810400515 0.0566663810400515 0.0566663810400515 0.773334475839794 0.0566663810400515]), 
    SVector{5}([0.0566663810400515 0.0566663810400515 0.773334475839794 0.0566663810400515 0.0566663810400515]), 
    SVector{5}([0.0566663810400515 0.773334475839794 0.0566663810400515 0.0566663810400515 0.0566663810400515]), 
    SVector{5}([0.773334475839794 0.0566663810400515 0.0566663810400515 0.0566663810400515 0.0566663810400515]), 
    SVector{5}([0.0828237846356080 0.0828237846356080 0.0828237846356080 0.375764323046588 0.375764323046588] ), 
    SVector{5}([0.0828237846356080 0.0828237846356080 0.375764323046588 0.0828237846356080 0.375764323046588] ), 
    SVector{5}([0.0828237846356080 0.0828237846356080 0.375764323046588 0.375764323046588 0.0828237846356080] ), 
    SVector{5}([0.0828237846356080 0.375764323046588 0.0828237846356080 0.0828237846356080 0.375764323046588] ), 
    SVector{5}([0.0828237846356080 0.375764323046588 0.0828237846356080 0.375764323046588 0.0828237846356080] ), 
    SVector{5}([0.0828237846356080 0.375764323046588 0.375764323046588 0.0828237846356080 0.0828237846356080] ), 
    SVector{5}([0.375764323046588 0.0828237846356080 0.0828237846356080 0.0828237846356080 0.375764323046588] ), 
    SVector{5}([0.375764323046588 0.0828237846356080 0.0828237846356080 0.375764323046588 0.0828237846356080] ), 
    SVector{5}([0.375764323046588 0.0828237846356080 0.375764323046588 0.0828237846356080 0.0828237846356080] ), 
    SVector{5}([0.375764323046588 0.375764323046588 0.0828237846356080 0.0828237846356080 0.0828237846356080] )
]

const pentatope15w = [
    0.019717445949776514,
    0.019717445949776514,
    0.019717445949776514,
    0.019717445949776514,
    0.019717445949776514,
    0.09014127702511174,
    0.09014127702511174,
    0.09014127702511174,
    0.09014127702511174,
    0.09014127702511174,
    0.09014127702511174,
    0.09014127702511174,
    0.09014127702511174,
    0.09014127702511174,
    0.09014127702511174
]


const pentatope35a = [

    SVector{5}([0.0863927292322510 0.0863927292322510 0.0863927292322510 0.0863927292322510 0.654429083070996]), 
    SVector{5}([0.0863927292322510 0.0863927292322510 0.0863927292322510 0.654429083070996 0.0863927292322510]), 
    SVector{5}([0.0863927292322510 0.0863927292322510 0.654429083070996 0.0863927292322510 0.0863927292322510]), 
    SVector{5}([0.0863927292322510 0.654429083070996 0.0863927292322510 0.0863927292322510 0.0863927292322510]), 
    SVector{5}([0.654429083070996 0.0863927292322510 0.0863927292322510 0.0863927292322510 0.0863927292322510]), 
    SVector{5}([0.0240149672006202 0.0240149672006202 0.0240149672006202 0.463977549199070 0.463977549199070] ), 
    SVector{5}([0.0240149672006202 0.0240149672006202 0.463977549199070 0.0240149672006202 0.463977549199070] ), 
    SVector{5}([0.0240149672006202 0.0240149672006202 0.463977549199070 0.463977549199070 0.0240149672006202] ), 
    SVector{5}([0.0240149672006202 0.463977549199070 0.0240149672006202 0.0240149672006202 0.463977549199070] ), 
    SVector{5}([0.0240149672006202 0.463977549199070 0.0240149672006202 0.463977549199070 0.0240149672006202] ), 
    SVector{5}([0.0240149672006202 0.463977549199070 0.463977549199070 0.0240149672006202 0.0240149672006202] ), 
    SVector{5}([0.463977549199070 0.0240149672006202 0.0240149672006202 0.0240149672006202 0.463977549199070] ), 
    SVector{5}([0.463977549199070 0.0240149672006202 0.0240149672006202 0.463977549199070 0.0240149672006202] ), 
    SVector{5}([0.463977549199070 0.0240149672006202 0.463977549199070 0.0240149672006202 0.0240149672006202] ), 
    SVector{5}([0.463977549199070 0.463977549199070 0.0240149672006202 0.0240149672006202 0.0240149672006202] ), 
    SVector{5}([0.293818004028937 0.293818004028937 0.293818004028937 0.0624751755625809 0.0560708123506085]  ), 
    SVector{5}([0.293818004028937 0.293818004028937 0.293818004028937 0.0560708123506085 0.0624751755625809]  ), 
    SVector{5}([0.293818004028937 0.293818004028937 0.0624751755625809 0.293818004028937 0.0560708123506085]  ), 
    SVector{5}([0.293818004028937 0.293818004028937 0.0624751755625809 0.0560708123506085 0.293818004028937]  ), 
    SVector{5}([0.293818004028937 0.293818004028937 0.0560708123506085 0.293818004028937 0.0624751755625809]  ), 
    SVector{5}([0.293818004028937 0.293818004028937 0.0560708123506085 0.0624751755625809 0.293818004028937]  ), 
    SVector{5}([0.293818004028937 0.0624751755625809 0.293818004028937 0.293818004028937 0.0560708123506085]  ), 
    SVector{5}([0.293818004028937 0.0624751755625809 0.293818004028937 0.0560708123506085 0.293818004028937]  ), 
    SVector{5}([0.293818004028937 0.0624751755625809 0.0560708123506085 0.293818004028937 0.293818004028937]  ), 
    SVector{5}([0.293818004028937 0.0560708123506085 0.293818004028937 0.293818004028937 0.0624751755625809]  ), 
    SVector{5}([0.293818004028937 0.0560708123506085 0.293818004028937 0.0624751755625809 0.293818004028937]  ), 
    SVector{5}([0.293818004028937 0.0560708123506085 0.0624751755625809 0.293818004028937 0.293818004028937]  ), 
    SVector{5}([0.0624751755625809 0.293818004028937 0.293818004028937 0.293818004028937 0.0560708123506085]  ), 
    SVector{5}([0.0624751755625809 0.293818004028937 0.293818004028937 0.0560708123506085 0.293818004028937]  ), 
    SVector{5}([0.0624751755625809 0.293818004028937 0.0560708123506085 0.293818004028937 0.293818004028937]  ), 
    SVector{5}([0.0624751755625809 0.0560708123506085 0.293818004028937 0.293818004028937 0.293818004028937]  ), 
    SVector{5}([0.0560708123506085 0.293818004028937 0.293818004028937 0.293818004028937 0.0624751755625809]  ), 
    SVector{5}([0.0560708123506085 0.293818004028937 0.293818004028937 0.0624751755625809 0.293818004028937]  ), 
    SVector{5}([0.0560708123506085 0.293818004028937 0.0624751755625809 0.293818004028937 0.293818004028937]  ), 
    SVector{5}([0.0560708123506085 0.0624751755625809 0.293818004028937 0.293818004028937 0.293818004028937]  )
]

const pentatope35w = [
    0.05144687284129604,
    0.05144687284129604,
    0.05144687284129604,
    0.05144687284129604,
    0.05144687284129604,
    0.010758106723188282,
    0.010758106723188282,
    0.010758106723188282,
    0.010758106723188282,
    0.010758106723188282,
    0.010758106723188282,
    0.010758106723188282,
    0.010758106723188282,
    0.010758106723188282,
    0.010758106723188282,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855,
    0.031759228428081855
]

const pentatopeA = [
  pentatope1a,
  pentatope5a,
  pentatope15a,
  pentatope35a
]

const pentatopeW = [
  pentatope1w,
  pentatope5w,
  pentatope15w,
  pentatope35w
]