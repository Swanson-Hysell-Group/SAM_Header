# Instructions for using SAM_Header to check block sample orientations

1. Our block samples are often collected with both magnetic compass orientations and sun compass orientations. However, the sun compass used for block sampling is different from the sun compass on a pomeroy aperture. The sun compass for block sampling has the tick marks going clockwise, whereas the pomeroy sun compass has the tick marks going counter clockwise.
2. We can still use the same SAM_Header Python script to check the consistency of block orientations collected in the field. The only change needed is during data recording - one needs to put in the negative number of what was recorded for the sun angle of the block!

