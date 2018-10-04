import numpy as np
import scipy

def truncNorm( a=0, b=1, mu=0, sigma=1, size=1 ):
    return scipy.stats.truncnorm.rvs( (a-mu)/sigma, (b-mu)/sigma, loc=mu, scale=sigma, size=size )

def truncNormCdf( x, a=0, b=1, mu=0, sigma=1 ):
    return scipy.stats.truncnorm.cdf( x, (a-mu)/sigma, (b-mu)/sigma, loc=mu, scale=sigma )

def truncNormPdf( x, a=0, b=1, mu=0, sigma=1 ):
    return scipy.stats.truncnorm.pdf( x, (a-mu)/sigma, (b-mu)/sigma, loc=mu, scale=sigma )

def truncPowerLawRVs( k, minX, maxX, size=1, skew='negative', algorithm='custom' ):
    if algorithm == 'custom':
        u = np.random.random( size )
        if skew == 'negative':
            vals = maxX - (maxX - minX)*np.power(1-u,1/(1-k))
        elif skew == 'positive':
            k = -k
            vals = np.power( (np.power(maxX,k+1) - np.power(minX,k+1))*u + np.power(minX,k+1), 1/(k+1) )
    elif algorithm == 'scipy':
        vals = scipy.stats.powerlaw.rvs( k+1, loc=minX, scale=(maxX-minX), size=size )  

    return vals

def truncPowerLawPdf( x, k, minX, maxX ):
    if not isinstance(x,np.ndarray):
        x = np.array( x )

    pdf = np.zeros_like( x )
    pdf[x < minX] = 0
    pdf[x > maxX] = 0
    remaining = (x > minX ) & ( x < maxX )
    pdf[remaining] = np.power(1/(maxX-x[remaining]),k)

    return pdf
    
def truncPowerLawMean( k, minX, maxX ):
    m = maxX - ((1-k)/(2-k))*(maxX - minX)
    return m

def truncPowerLawVar( k, minX, maxX ):
    eDiff = maxX - minX
    secondmoment = (1 - k) * ( eDiff*eDiff/(3-k) + maxX*maxX/(1-k) - 2*maxX*eDiff/(2-k) )
    #eDiff = maxX
    #print( (1- k) * ( eDiff*eDiff/(3-k) + maxX*maxX/(1-k) - 2*maxX*eDiff/(2-k) )- truncPowerLawMean( k, 0, maxX )**2  )
    #print( maxX * maxX * (1 + (k-4)*(1-k)/((3-k)*(2-k))) - truncPowerLawMean( k, 0, maxX )**2 )
    return secondmoment - truncPowerLawMean( k, minX, maxX )**2


