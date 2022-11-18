\name{generationsDistance}
\alias{generationsDistance}
\alias{generationsDistance-methods}
\alias{generationsDistance,ANY-method}

\title{ Calculate the distance between 2 generations of cells on the FACS scale }

\description{
  This function calculate the distance between 2 generations of cells on the FACS scale.
}

\usage{
  generationsDistance(dataRange, logDecades)
}

\arguments{
  \item{dataRange}{
    Digital Data range on the FACS instrument
  }
  \item{logDecades}{
    Number of log decades on the FACS instrument (dynamic range)
  }
}

\value{
  Return the spacing between generations on the FACS scale.
}

\details{
  We can use this formula to convert FFI (FACS fluorescence Intensity) to RFI (Relative Fluorescence Intensity): \eqn{ RFI = 10^{\frac {FFI \cdot l}{c}} }{ RFI = 10^((FFI * l)/c) }

  The inverse formula is used to convert from RFI to FACS fluorescence: \eqn{ FFI = \frac{c*\log(RFI)}{(l*\log(10))} }{ FFI = (c*log(RFI)) / (l*log(10)) }

  Where:
  \describe{
    \item{RFI}{is the Relative Fluorescence Intensity}
    \item{FFI}{is the fluorescence on the FACS scale}
    \item{\eqn{l}}{is the number of log decades in the FACS instrument}
    \item{\eqn{c}}{is the number of data points (channels) in the instrument.}
  }

  Using this formulas it is possible to estimate the spacing between generations on the FACS scale.
  The spacing value is automatically computed, based on the number of decades and the assumption that each generation has one-half of the intensity of the previous generation.
}

\author{ Davide Rambaldi}

\examples{
  distance <- generationsDistance(1024, 4)
}

\references{
  \enumerate{
    \item \bold{Tracking antigen-driven responses by flow cytometry: Monitoring proliferation by dye dilution.} Paul K Wallace, Joseph D Tario, Jan L Fisher, Stephen S Wallace, Marc S Ernstoff, Katharine A Muirhead. Cytometry (2008) vol. 73A (11) pp. 1019-1034
    \item FACS Formulas to convert between Relative Fluorescence and fluorescence on the FACS scale: \bold{MEASURING MOLECULES OF EQUIVALENT FLUOROCHROME (MEF) USING SPHEROTM RAINBOW AND ULTRA RAINBOW CALIBRATION PARTICLES}, Spherotech, http://www.spherotech.com/tech_SpheroTech_Note_9.html
  }
}