!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_sun
!
! Sun zenith and declination
!
   use mod_intkinds
   use mod_realkinds
   use mod_stdio
   use mod_memutil
   use mod_constants
   use mod_dynparam
   use mod_runparams
   use mod_mpmessage
   use mod_mppparam
   use mod_service
   use mod_date
   use mod_sunorbit
   use netcdf

   implicit none

   private

   real(rkx) , dimension(:) , pointer :: heppatsi => null()
   real(rkx) , dimension(3,1610:2008) :: tsi
   real(rkx) , parameter :: tsifac = 0.9965_rkx
   integer(ik4) :: ii , jj

   public :: zenitm
   !
   ! ANNUAL MEAN TSI
   ! Lean (GRL 2000) with Wang Lean Sheeley (ApJ 2005) background
   ! Mon Apr  6 11:29:27 2009
   ! PMOD absolute scale , multiply by 0.9965 for TIM scale
   ! YEAR, TSI 11-yr cycle, TSI cycle+background
   !
   data ((tsi(ii,jj),ii=1,3),jj=1610,2008) / &
       1610.5_rkx,1365.8477_rkx,1365.5469_rkx, 1611.5_rkx,1365.8342_rkx,1365.5300_rkx, &
       1612.5_rkx,1366.2461_rkx,1365.9279_rkx, 1613.5_rkx,1366.3650_rkx,1366.0399_rkx, &
       1614.5_rkx,1366.4451_rkx,1366.1143_rkx, 1615.5_rkx,1366.1591_rkx,1365.8314_rkx, &
       1616.5_rkx,1365.7358_rkx,1365.4148_rkx, 1617.5_rkx,1365.6107_rkx,1365.2889_rkx, &
       1618.5_rkx,1365.6038_rkx,1365.2783_rkx, 1619.5_rkx,1365.7001_rkx,1365.3684_rkx, &
       1620.5_rkx,1365.7001_rkx,1365.3645_rkx, 1621.5_rkx,1365.7001_rkx,1365.3607_rkx, &
       1622.5_rkx,1365.7001_rkx,1365.3568_rkx, 1623.5_rkx,1365.7001_rkx,1365.3530_rkx, &
       1624.5_rkx,1365.6621_rkx,1365.3121_rkx, 1625.5_rkx,1365.8926_rkx,1365.5303_rkx, &
       1626.5_rkx,1365.7816_rkx,1365.4191_rkx, 1627.5_rkx,1365.7106_rkx,1365.3418_rkx, &
       1628.5_rkx,1365.7577_rkx,1365.3518_rkx, 1629.5_rkx,1365.7261_rkx,1365.2922_rkx, &
       1630.5_rkx,1365.5946_rkx,1365.1428_rkx, 1631.5_rkx,1365.6255_rkx,1365.1515_rkx, &
       1632.5_rkx,1365.5946_rkx,1365.1183_rkx, 1633.5_rkx,1365.6951_rkx,1365.2158_rkx, &
       1634.5_rkx,1365.6157_rkx,1365.1362_rkx, 1635.5_rkx,1365.6249_rkx,1365.1411_rkx, &
       1636.5_rkx,1365.5946_rkx,1365.1080_rkx, 1637.5_rkx,1365.5946_rkx,1365.1046_rkx, &
       1638.5_rkx,1366.0768_rkx,1365.5710_rkx, 1639.5_rkx,1366.1344_rkx,1365.6241_rkx, &
       1640.5_rkx,1365.7001_rkx,1365.1936_rkx, 1641.5_rkx,1365.5946_rkx,1365.0815_rkx, &
       1642.5_rkx,1365.9277_rkx,1365.4006_rkx, 1643.5_rkx,1365.7183_rkx,1365.1824_rkx, &
       1644.5_rkx,1365.6761_rkx,1365.1272_rkx, 1645.5_rkx,1365.5946_rkx,1365.0454_rkx, &
       1646.5_rkx,1365.5946_rkx,1365.0449_rkx, 1647.5_rkx,1365.5946_rkx,1365.0443_rkx, &
       1648.5_rkx,1365.5946_rkx,1365.0424_rkx, 1649.5_rkx,1365.5946_rkx,1365.0399_rkx, &
       1650.5_rkx,1365.5946_rkx,1365.0389_rkx, 1651.5_rkx,1365.5946_rkx,1365.0383_rkx, &
       1652.5_rkx,1365.6227_rkx,1365.0657_rkx, 1653.5_rkx,1365.6010_rkx,1365.0439_rkx, &
       1654.5_rkx,1365.5995_rkx,1365.0358_rkx, 1655.5_rkx,1365.5981_rkx,1365.0260_rkx, &
       1656.5_rkx,1365.5989_rkx,1365.0249_rkx, 1657.5_rkx,1365.5961_rkx,1365.0199_rkx, &
       1658.5_rkx,1365.5946_rkx,1365.0145_rkx, 1659.5_rkx,1365.5946_rkx,1365.0125_rkx, &
       1660.5_rkx,1365.6086_rkx,1365.0259_rkx, 1661.5_rkx,1365.6002_rkx,1365.0178_rkx, &
       1662.5_rkx,1365.5946_rkx,1365.0125_rkx, 1663.5_rkx,1365.5946_rkx,1365.0125_rkx, &
       1664.5_rkx,1365.5946_rkx,1365.0126_rkx, 1665.5_rkx,1365.5946_rkx,1365.0127_rkx, &
       1666.5_rkx,1365.5946_rkx,1365.0127_rkx, 1667.5_rkx,1365.5946_rkx,1365.0125_rkx, &
       1668.5_rkx,1365.5946_rkx,1365.0122_rkx, 1669.5_rkx,1365.5946_rkx,1365.0122_rkx, &
       1670.5_rkx,1365.5946_rkx,1365.0122_rkx, 1671.5_rkx,1365.6010_rkx,1365.0183_rkx, &
       1672.5_rkx,1365.5974_rkx,1365.0148_rkx, 1673.5_rkx,1365.5946_rkx,1365.0122_rkx, &
       1674.5_rkx,1365.5961_rkx,1365.0135_rkx, 1675.5_rkx,1365.5946_rkx,1365.0120_rkx, &
       1676.5_rkx,1365.6073_rkx,1365.0239_rkx, 1677.5_rkx,1365.5967_rkx,1365.0134_rkx, &
       1678.5_rkx,1365.5961_rkx,1365.0128_rkx, 1679.5_rkx,1365.5946_rkx,1365.0115_rkx, &
       1680.5_rkx,1365.6002_rkx,1365.0170_rkx, 1681.5_rkx,1365.5946_rkx,1365.0115_rkx, &
       1682.5_rkx,1365.5946_rkx,1365.0115_rkx, 1683.5_rkx,1365.5946_rkx,1365.0115_rkx, &
       1684.5_rkx,1365.6045_rkx,1365.0211_rkx, 1685.5_rkx,1365.5946_rkx,1365.0117_rkx, &
       1686.5_rkx,1365.5989_rkx,1365.0159_rkx, 1687.5_rkx,1365.5953_rkx,1365.0126_rkx, &
       1688.5_rkx,1365.5981_rkx,1365.0159_rkx, 1689.5_rkx,1365.5961_rkx,1365.0146_rkx, &
       1690.5_rkx,1365.5946_rkx,1365.0143_rkx, 1691.5_rkx,1365.5946_rkx,1365.0157_rkx, &
       1692.5_rkx,1365.5946_rkx,1365.0172_rkx, 1693.5_rkx,1365.5946_rkx,1365.0177_rkx, &
       1694.5_rkx,1365.5946_rkx,1365.0179_rkx, 1695.5_rkx,1365.5953_rkx,1365.0186_rkx, &
       1696.5_rkx,1365.5946_rkx,1365.0178_rkx, 1697.5_rkx,1365.5946_rkx,1365.0178_rkx, &
       1698.5_rkx,1365.5946_rkx,1365.0179_rkx, 1699.5_rkx,1365.5946_rkx,1365.0184_rkx, &
       1700.5_rkx,1365.5974_rkx,1365.0216_rkx, 1701.5_rkx,1365.5981_rkx,1365.0236_rkx, &
       1702.5_rkx,1365.5989_rkx,1365.0266_rkx, 1703.5_rkx,1365.6135_rkx,1365.0444_rkx, &
       1704.5_rkx,1365.6234_rkx,1365.0594_rkx, 1705.5_rkx,1365.6333_rkx,1365.0752_rkx, &
       1706.5_rkx,1365.6171_rkx,1365.0637_rkx, 1707.5_rkx,1365.6318_rkx,1365.0802_rkx, &
       1708.5_rkx,1365.6143_rkx,1365.0658_rkx, 1709.5_rkx,1365.6058_rkx,1365.0614_rkx, &
       1710.5_rkx,1365.5974_rkx,1365.0634_rkx, 1711.5_rkx,1365.5946_rkx,1365.0739_rkx, &
       1712.5_rkx,1365.5946_rkx,1365.0798_rkx, 1713.5_rkx,1365.5967_rkx,1365.0863_rkx, &
       1714.5_rkx,1365.6010_rkx,1365.1023_rkx, 1715.5_rkx,1365.6199_rkx,1365.1294_rkx, &
       1716.5_rkx,1365.6586_rkx,1365.1694_rkx, 1717.5_rkx,1365.7177_rkx,1365.2294_rkx, &
       1718.5_rkx,1365.6578_rkx,1365.1707_rkx, 1719.5_rkx,1365.8329_rkx,1365.3429_rkx, &
       1720.5_rkx,1365.7590_rkx,1365.2859_rkx, 1721.5_rkx,1365.7450_rkx,1365.2880_rkx, &
       1722.5_rkx,1365.6719_rkx,1365.2209_rkx, 1723.5_rkx,1365.6255_rkx,1365.1837_rkx, &
       1724.5_rkx,1365.7042_rkx,1365.2681_rkx, 1725.5_rkx,1365.6846_rkx,1365.2574_rkx, &
       1726.5_rkx,1365.8490_rkx,1365.4274_rkx, 1727.5_rkx,1365.8512_rkx,1365.4327_rkx, &
       1728.5_rkx,1366.0459_rkx,1365.6237_rkx, 1729.5_rkx,1365.7633_rkx,1365.3479_rkx, &
       1730.5_rkx,1366.0845_rkx,1365.6605_rkx, 1731.5_rkx,1365.5946_rkx,1365.1812_rkx, &
       1732.5_rkx,1365.7211_rkx,1365.3090_rkx, 1733.5_rkx,1365.5946_rkx,1365.1984_rkx, &
       1734.5_rkx,1365.5946_rkx,1365.2086_rkx, 1735.5_rkx,1365.7233_rkx,1365.3386_rkx, &
       1736.5_rkx,1365.9362_rkx,1365.5486_rkx, 1737.5_rkx,1365.7633_rkx,1365.3827_rkx, &
       1738.5_rkx,1365.7141_rkx,1365.3370_rkx, 1739.5_rkx,1365.9636_rkx,1365.5795_rkx, &
       1740.5_rkx,1365.6599_rkx,1365.2811_rkx, 1741.5_rkx,1366.0001_rkx,1365.6107_rkx, &
       1742.5_rkx,1365.7078_rkx,1365.3247_rkx, 1743.5_rkx,1365.6530_rkx,1365.2698_rkx, &
       1744.5_rkx,1365.5946_rkx,1365.2137_rkx, 1745.5_rkx,1365.5946_rkx,1365.2164_rkx, &
       1746.5_rkx,1365.5946_rkx,1365.2223_rkx, 1747.5_rkx,1365.5946_rkx,1365.2305_rkx, &
       1748.5_rkx,1366.0233_rkx,1365.6548_rkx, 1749.5_rkx,1366.0388_rkx,1365.6749_rkx, &
       1750.5_rkx,1366.0023_rkx,1365.6385_rkx, 1751.5_rkx,1365.8314_rkx,1365.4680_rkx, &
       1752.5_rkx,1365.7985_rkx,1365.4402_rkx, 1753.5_rkx,1365.7626_rkx,1365.4220_rkx, &
       1754.5_rkx,1365.6565_rkx,1365.3351_rkx, 1755.5_rkx,1365.6277_rkx,1365.3220_rkx, &
       1756.5_rkx,1365.6459_rkx,1365.3501_rkx, 1757.5_rkx,1365.7689_rkx,1365.4734_rkx, &
       1758.5_rkx,1365.8806_rkx,1365.5867_rkx, 1759.5_rkx,1365.9425_rkx,1365.6501_rkx, &
       1760.5_rkx,1365.9144_rkx,1365.6254_rkx, 1761.5_rkx,1366.0760_rkx,1365.7898_rkx, &
       1762.5_rkx,1365.9193_rkx,1365.6514_rkx, 1763.5_rkx,1365.8350_rkx,1365.5811_rkx, &
       1764.5_rkx,1365.8090_rkx,1365.5573_rkx, 1765.5_rkx,1365.6537_rkx,1365.4065_rkx, &
       1766.5_rkx,1365.6206_rkx,1365.3759_rkx, 1767.5_rkx,1365.8329_rkx,1365.5817_rkx, &
       1768.5_rkx,1366.0957_rkx,1365.8346_rkx, 1769.5_rkx,1366.2869_rkx,1366.0194_rkx, &
       1770.5_rkx,1366.2806_rkx,1366.0220_rkx, 1771.5_rkx,1366.1527_rkx,1365.9155_rkx, &
       1772.5_rkx,1366.0599_rkx,1365.8433_rkx, 1773.5_rkx,1365.8224_rkx,1365.6243_rkx, &
       1774.5_rkx,1365.7760_rkx,1365.5862_rkx, 1775.5_rkx,1365.6339_rkx,1365.4493_rkx, &
       1776.5_rkx,1365.6937_rkx,1365.5039_rkx, 1777.5_rkx,1365.8638_rkx,1365.6656_rkx, &
       1778.5_rkx,1366.1007_rkx,1365.8955_rkx, 1779.5_rkx,1366.1625_rkx,1365.9534_rkx, &
       1780.5_rkx,1365.9812_rkx,1365.7753_rkx, 1781.5_rkx,1366.0944_rkx,1365.8868_rkx, &
       1782.5_rkx,1365.8258_rkx,1365.6217_rkx, 1783.5_rkx,1365.7429_rkx,1365.5284_rkx, &
       1784.5_rkx,1365.6283_rkx,1365.3966_rkx, 1785.5_rkx,1365.7070_rkx,1365.4558_rkx, &
       1786.5_rkx,1366.0396_rkx,1365.7682_rkx, 1787.5_rkx,1366.2216_rkx,1365.9337_rkx, &
       1788.5_rkx,1366.1744_rkx,1365.8801_rkx, 1789.5_rkx,1366.1548_rkx,1365.8595_rkx, &
       1790.5_rkx,1366.0521_rkx,1365.7604_rkx, 1791.5_rkx,1365.8982_rkx,1365.6101_rkx, &
       1792.5_rkx,1365.8898_rkx,1365.5961_rkx, 1793.5_rkx,1365.8828_rkx,1365.5756_rkx, &
       1794.5_rkx,1365.8069_rkx,1365.4816_rkx, 1795.5_rkx,1365.7050_rkx,1365.3645_rkx, &
       1796.5_rkx,1365.6909_rkx,1365.3348_rkx, 1797.5_rkx,1365.6487_rkx,1365.2817_rkx, &
       1798.5_rkx,1365.6277_rkx,1365.2567_rkx, 1799.5_rkx,1365.6339_rkx,1365.2629_rkx, &
       1800.5_rkx,1365.6719_rkx,1365.3035_rkx, 1801.5_rkx,1365.9537_rkx,1365.5757_rkx, &
       1802.5_rkx,1365.8428_rkx,1365.4541_rkx, 1803.5_rkx,1365.7246_rkx,1365.3218_rkx, &
       1804.5_rkx,1365.7465_rkx,1365.3257_rkx, 1805.5_rkx,1365.7745_rkx,1365.3361_rkx, &
       1806.5_rkx,1365.6881_rkx,1365.2385_rkx, 1807.5_rkx,1365.6298_rkx,1365.1710_rkx, &
       1808.5_rkx,1365.6193_rkx,1365.1458_rkx, 1809.5_rkx,1365.6030_rkx,1365.1180_rkx, &
       1810.5_rkx,1365.5946_rkx,1365.1094_rkx, 1811.5_rkx,1365.5967_rkx,1365.1222_rkx, &
       1812.5_rkx,1365.6227_rkx,1365.1631_rkx, 1813.5_rkx,1365.6586_rkx,1365.2117_rkx, &
       1814.5_rkx,1365.6677_rkx,1365.2355_rkx, 1815.5_rkx,1365.7126_rkx,1365.2906_rkx, &
       1816.5_rkx,1365.8110_rkx,1365.3866_rkx, 1817.5_rkx,1365.7914_rkx,1365.3600_rkx, &
       1818.5_rkx,1365.7471_rkx,1365.3119_rkx, 1819.5_rkx,1365.7295_rkx,1365.2968_rkx, &
       1820.5_rkx,1365.6698_rkx,1365.2410_rkx, 1821.5_rkx,1365.6249_rkx,1365.2194_rkx, &
       1822.5_rkx,1365.6157_rkx,1365.2432_rkx, 1823.5_rkx,1365.6030_rkx,1365.2483_rkx, &
       1824.5_rkx,1365.6305_rkx,1365.2893_rkx, 1825.5_rkx,1365.6958_rkx,1365.3627_rkx, &
       1826.5_rkx,1365.7957_rkx,1365.4659_rkx, 1827.5_rkx,1365.9066_rkx,1365.5771_rkx, &
       1828.5_rkx,1365.9952_rkx,1365.6646_rkx, 1829.5_rkx,1366.0107_rkx,1365.6825_rkx, &
       1830.5_rkx,1366.0465_rkx,1365.7235_rkx, 1831.5_rkx,1365.8701_rkx,1365.5585_rkx, &
       1832.5_rkx,1365.7542_rkx,1365.4565_rkx, 1833.5_rkx,1365.6403_rkx,1365.3612_rkx, &
       1834.5_rkx,1365.6635_rkx,1365.3966_rkx, 1835.5_rkx,1365.9200_rkx,1365.6577_rkx, &
       1836.5_rkx,1366.2988_rkx,1366.0393_rkx, 1837.5_rkx,1366.3635_rkx,1366.1127_rkx, &
       1838.5_rkx,1366.1387_rkx,1365.8997_rkx, 1839.5_rkx,1366.0521_rkx,1365.8174_rkx, &
       1840.5_rkx,1365.9341_rkx,1365.7007_rkx, 1841.5_rkx,1365.7816_rkx,1365.5490_rkx, &
       1842.5_rkx,1365.7267_rkx,1365.4940_rkx, 1843.5_rkx,1365.6522_rkx,1365.4237_rkx, &
       1844.5_rkx,1365.6775_rkx,1365.4543_rkx, 1845.5_rkx,1365.8041_rkx,1365.5820_rkx, &
       1846.5_rkx,1365.9025_rkx,1365.6803_rkx, 1847.5_rkx,1366.0023_rkx,1365.7816_rkx, &
       1848.5_rkx,1366.1976_rkx,1365.9781_rkx, 1849.5_rkx,1366.1829_rkx,1365.9691_rkx, &
       1850.5_rkx,1365.9812_rkx,1365.7698_rkx, 1851.5_rkx,1366.0029_rkx,1365.7744_rkx, &
       1852.5_rkx,1365.9446_rkx,1365.6947_rkx, 1853.5_rkx,1365.8448_rkx,1365.5848_rkx, &
       1854.5_rkx,1365.7162_rkx,1365.4614_rkx, 1855.5_rkx,1365.6262_rkx,1365.3828_rkx, &
       1856.5_rkx,1365.6163_rkx,1365.3853_rkx, 1857.5_rkx,1365.7169_rkx,1365.4946_rkx, &
       1858.5_rkx,1365.9066_rkx,1365.6875_rkx, 1859.5_rkx,1366.1260_rkx,1365.9054_rkx, &
       1860.5_rkx,1366.1963_rkx,1365.9718_rkx, 1861.5_rkx,1366.0916_rkx,1365.8623_rkx, &
       1862.5_rkx,1365.9496_rkx,1365.7120_rkx, 1863.5_rkx,1365.8821_rkx,1365.6283_rkx, &
       1864.5_rkx,1365.8370_rkx,1365.5660_rkx, 1865.5_rkx,1365.7534_rkx,1365.4755_rkx, &
       1866.5_rkx,1365.6909_rkx,1365.4119_rkx, 1867.5_rkx,1365.6382_rkx,1365.3597_rkx, &
       1868.5_rkx,1365.7977_rkx,1365.5194_rkx, 1869.5_rkx,1366.0325_rkx,1365.7557_rkx, &
       1870.5_rkx,1366.2708_rkx,1365.9944_rkx, 1871.5_rkx,1366.2054_rkx,1365.9343_rkx, &
       1872.5_rkx,1366.1576_rkx,1365.8876_rkx, 1873.5_rkx,1365.9580_rkx,1365.6866_rkx, &
       1874.5_rkx,1365.8406_rkx,1365.5582_rkx, 1875.5_rkx,1365.7035_rkx,1365.4095_rkx, &
       1876.5_rkx,1365.6586_rkx,1365.3593_rkx, 1877.5_rkx,1365.6543_rkx,1365.3596_rkx, &
       1878.5_rkx,1365.6135_rkx,1365.3309_rkx, 1879.5_rkx,1365.6255_rkx,1365.3533_rkx, &
       1880.5_rkx,1365.7689_rkx,1365.5000_rkx, 1881.5_rkx,1365.9124_rkx,1365.6443_rkx, &
       1882.5_rkx,1365.9313_rkx,1365.6676_rkx, 1883.5_rkx,1365.9791_rkx,1365.7147_rkx, &
       1884.5_rkx,1365.8812_rkx,1365.6166_rkx, 1885.5_rkx,1365.7909_rkx,1365.5070_rkx, &
       1886.5_rkx,1365.6487_rkx,1365.3417_rkx, 1887.5_rkx,1365.6234_rkx,1365.2982_rkx, &
       1888.5_rkx,1365.5962_rkx,1365.2628_rkx, 1889.5_rkx,1365.5652_rkx,1365.2344_rkx, &
       1890.5_rkx,1365.5912_rkx,1365.2690_rkx, 1891.5_rkx,1365.8303_rkx,1365.5204_rkx, &
       1892.5_rkx,1365.9163_rkx,1365.6190_rkx, 1893.5_rkx,1366.0458_rkx,1365.7600_rkx, &
       1894.5_rkx,1366.1332_rkx,1365.8553_rkx, 1895.5_rkx,1366.0166_rkx,1365.7390_rkx, &
       1896.5_rkx,1365.8434_rkx,1365.5581_rkx, 1897.5_rkx,1365.7094_rkx,1365.4126_rkx, &
       1898.5_rkx,1365.6982_rkx,1365.3899_rkx, 1899.5_rkx,1365.6534_rkx,1365.3381_rkx, &
       1900.5_rkx,1365.6216_rkx,1365.3074_rkx, 1901.5_rkx,1365.5294_rkx,1365.2292_rkx, &
       1902.5_rkx,1365.5165_rkx,1365.2378_rkx, 1903.5_rkx,1365.7083_rkx,1365.4479_rkx, &
       1904.5_rkx,1365.9651_rkx,1365.7180_rkx, 1905.5_rkx,1365.7684_rkx,1365.5291_rkx, &
       1906.5_rkx,1365.9651_rkx,1365.7255_rkx, 1907.5_rkx,1365.8604_rkx,1365.6097_rkx, &
       1908.5_rkx,1365.9426_rkx,1365.6748_rkx, 1909.5_rkx,1365.8459_rkx,1365.5642_rkx, &
       1910.5_rkx,1365.7173_rkx,1365.4309_rkx, 1911.5_rkx,1365.6285_rkx,1365.3473_rkx, &
       1912.5_rkx,1365.5706_rkx,1365.3010_rkx, 1913.5_rkx,1365.5739_rkx,1365.3175_rkx, &
       1914.5_rkx,1365.6302_rkx,1365.3844_rkx, 1915.5_rkx,1365.9285_rkx,1365.6890_rkx, &
       1916.5_rkx,1366.1349_rkx,1365.8990_rkx, 1917.5_rkx,1366.2821_rkx,1366.0480_rkx, &
       1918.5_rkx,1366.2454_rkx,1366.0096_rkx, 1919.5_rkx,1366.0179_rkx,1365.7802_rkx, &
       1920.5_rkx,1365.8523_rkx,1365.6178_rkx, 1921.5_rkx,1365.7351_rkx,1365.5127_rkx, &
       1922.5_rkx,1365.6019_rkx,1365.3948_rkx, 1923.5_rkx,1365.6211_rkx,1365.4265_rkx, &
       1924.5_rkx,1365.6436_rkx,1365.4581_rkx, 1925.5_rkx,1365.8406_rkx,1365.6622_rkx, &
       1926.5_rkx,1365.9348_rkx,1365.7633_rkx, 1927.5_rkx,1366.1135_rkx,1365.9468_rkx, &
       1928.5_rkx,1365.9885_rkx,1365.8245_rkx, 1929.5_rkx,1365.9429_rkx,1365.7833_rkx, &
       1930.5_rkx,1365.9159_rkx,1365.7655_rkx, 1931.5_rkx,1365.7780_rkx,1365.6436_rkx, &
       1932.5_rkx,1365.6583_rkx,1365.5364_rkx, 1933.5_rkx,1365.5300_rkx,1365.4156_rkx, &
       1934.5_rkx,1365.6361_rkx,1365.5275_rkx, 1935.5_rkx,1365.8500_rkx,1365.7439_rkx, &
       1936.5_rkx,1366.2373_rkx,1366.1333_rkx, 1937.5_rkx,1366.1718_rkx,1366.0676_rkx, &
       1938.5_rkx,1366.1079_rkx,1366.0031_rkx, 1939.5_rkx,1366.0894_rkx,1365.9868_rkx, &
       1940.5_rkx,1366.0143_rkx,1365.9242_rkx, 1941.5_rkx,1365.9130_rkx,1365.8451_rkx, &
       1942.5_rkx,1365.7847_rkx,1365.7419_rkx, 1943.5_rkx,1365.6052_rkx,1365.5841_rkx, &
       1944.5_rkx,1365.6224_rkx,1365.6140_rkx, 1945.5_rkx,1365.8850_rkx,1365.8810_rkx, &
       1946.5_rkx,1365.9818_rkx,1365.9791_rkx, 1947.5_rkx,1366.2190_rkx,1366.2185_rkx, &
       1948.5_rkx,1366.3475_rkx,1366.3490_rkx, 1949.5_rkx,1366.2528_rkx,1366.2555_rkx, &
       1950.5_rkx,1366.0098_rkx,1366.0131_rkx, 1951.5_rkx,1365.7721_rkx,1365.7765_rkx, &
       1952.5_rkx,1365.7653_rkx,1365.7676_rkx, 1953.5_rkx,1365.6313_rkx,1365.6284_rkx, &
       1954.5_rkx,1365.6599_rkx,1365.6564_rkx, 1955.5_rkx,1365.7793_rkx,1365.7773_rkx, &
       1956.5_rkx,1366.3097_rkx,1366.3109_rkx, 1957.5_rkx,1366.6632_rkx,1366.6681_rkx, &
       1958.5_rkx,1366.6246_rkx,1366.6328_rkx, 1959.5_rkx,1366.3717_rkx,1366.3828_rkx, &
       1960.5_rkx,1366.2682_rkx,1366.2767_rkx, 1961.5_rkx,1365.9230_rkx,1365.9199_rkx, &
       1962.5_rkx,1365.7656_rkx,1365.7484_rkx, 1963.5_rkx,1365.7152_rkx,1365.6963_rkx, &
       1964.5_rkx,1365.7114_rkx,1365.6976_rkx, 1965.5_rkx,1365.7378_rkx,1365.7341_rkx, &
       1966.5_rkx,1365.9058_rkx,1365.9178_rkx, 1967.5_rkx,1366.0889_rkx,1366.1143_rkx, &
       1968.5_rkx,1366.1295_rkx,1366.1644_rkx, 1969.5_rkx,1366.2069_rkx,1366.2476_rkx, &
       1970.5_rkx,1366.2036_rkx,1366.2426_rkx, 1971.5_rkx,1365.9354_rkx,1365.9580_rkx, &
       1972.5_rkx,1366.0519_rkx,1366.0525_rkx, 1973.5_rkx,1365.8131_rkx,1365.7991_rkx, &
       1974.5_rkx,1365.7448_rkx,1365.7271_rkx, 1975.5_rkx,1365.5466_rkx,1365.5345_rkx, &
       1976.5_rkx,1365.6458_rkx,1365.6453_rkx, 1977.5_rkx,1365.8248_rkx,1365.8331_rkx, &
       1978.5_rkx,1366.2616_rkx,1366.2747_rkx, 1979.5_rkx,1366.6193_rkx,1366.6348_rkx, &
       1980.5_rkx,1366.6323_rkx,1366.6482_rkx, 1981.5_rkx,1366.6829_rkx,1366.6951_rkx, &
       1982.5_rkx,1366.2808_rkx,1366.2859_rkx, 1983.5_rkx,1366.1989_rkx,1366.1992_rkx, &
       1984.5_rkx,1365.8088_rkx,1365.8103_rkx, 1985.5_rkx,1365.6382_rkx,1365.6416_rkx, &
       1986.5_rkx,1365.6345_rkx,1365.6379_rkx, 1987.5_rkx,1365.7865_rkx,1365.7899_rkx, &
       1988.5_rkx,1366.0792_rkx,1366.0826_rkx, 1989.5_rkx,1366.6445_rkx,1366.6479_rkx, &
       1990.5_rkx,1366.5499_rkx,1366.5533_rkx, 1991.5_rkx,1366.4423_rkx,1366.4457_rkx, &
       1992.5_rkx,1366.2987_rkx,1366.3021_rkx, 1993.5_rkx,1366.0251_rkx,1366.0286_rkx, &
       1994.5_rkx,1365.7937_rkx,1365.7971_rkx, 1995.5_rkx,1365.6962_rkx,1365.6996_rkx, &
       1996.5_rkx,1365.6086_rkx,1365.6121_rkx, 1997.5_rkx,1365.7365_rkx,1365.7399_rkx, &
       1998.5_rkx,1366.0986_rkx,1366.1021_rkx, 1999.5_rkx,1366.3817_rkx,1366.3851_rkx, &
       2000.5_rkx,1366.6836_rkx,1366.6836_rkx, 2001.5_rkx,1366.6022_rkx,1366.6022_rkx, &
       2002.5_rkx,1366.6807_rkx,1366.6807_rkx, 2003.5_rkx,1366.2300_rkx,1366.2300_rkx, &
       2004.5_rkx,1366.0480_rkx,1366.0480_rkx, 2005.5_rkx,1365.8545_rkx,1365.8545_rkx, &
       2006.5_rkx,1365.8107_rkx,1365.8107_rkx, 2007.5_rkx,1365.7240_rkx,1365.7240_rkx, &
       2008.5_rkx,1365.6918_rkx,1365.6918 /

  contains
  !
  ! This subroutine computes the solar declination angle
  ! from the julian date.
  !
  subroutine solar1
    implicit none
    real(rk8) :: decdeg , obliq , mvelp
    integer(ik4) , save :: lyear = bigint
    integer(ik4) :: iyear
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'solar1'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    iyear = rcmtimer%year + int(year_offset,ik4)
    if ( lyear /= iyear ) then
      if ( myid == italk ) then
        write (stdout, *) 'Updating orbit parameters at ',rcmtimer%str( )
      end if
      call orb_params(iyear,eccen,obliq,mvelp,obliqr,lambm0,mvelpp)
      lyear = iyear
    end if
    call orb_decl(yearpoint(rcmtimer%idate),eccen,mvelpp,lambm0, &
                  obliqr,declin,eccf)
    ! If we are fixing the solar constant, then fix the declination
    ! angle and eccentricity factor to constant values
    if ( irceideal == 1 ) then
      declin = 0.0_rkx
      eccf = 1.0_rkx
    end if
    decdeg = declin/degrad
    if ( myid == italk .and. alarm_day%act( ) ) then
      write (stdout, *) 'At ',rcmtimer%str( )
      write (stdout,'(a,f12.5,a,f12.8,a)') ' JDay ', calday , &
        ' solar declination angle = ', decdeg , ' degrees'
      write(stdout, '(18x,a,f12.4,a)') ' solar TSI irradiance    = ' , &
        solcon, ' W/m^2'
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine solar1
  !
  ! This subroutine calculates the cosine of the solar zenith angle
  ! for all longitude points of the RegCM domain. It needs as inputs
  ! the longitude and latitude of the points, the initial date of the
  ! simulation and the gmt. All these quantities are specified
  ! in the initialization procedure of RegCM
  !
  subroutine zenitm(xlat,xlon,coszrs)
    implicit none
    real(rkx) , pointer , intent(in), dimension(:,:) :: xlat , xlon
    real(rkx) , pointer , intent(inout), dimension(:,:) :: coszrs
    integer(ik4) :: i , j
    real(rkx) :: xxlat , xxlon
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'zenitm'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! Update solar constant for today
    !
    calday = real(yeardayfrac(rcmtimer%idate),rkx)
    if ( rcmtimer%start( ) .or. alarm_day%act( ) .or. doing_restart ) then
      if ( ifixsolar == 1 ) then
        ! Fix the solar constant; no diurnal or seasonal variability
        solcon = fixedsolarval
      else
        solcon = solar_irradiance( )
      end if
      call solar1( )
    end if
    scon = solcon*d_1000
    if ( ifixsolar == 1 ) then
      coszrs(:,:) = 1.0
    else
      do i = ici1 , ici2
        do j = jci1 , jci2
          xxlat = xlat(j,i)*degrad
          xxlon = xlon(j,i)*degrad
          coszrs(j,i) = orb_cosz(calday,xxlat,xxlon,declin)
          coszrs(j,i) = max(0.0_rkx,coszrs(j,i))
          coszrs(j,i) = min(1.0_rkx,coszrs(j,i))
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine zenitm

  real(rkx) function solar_irradiance( )
    implicit none
    integer(ik4) :: iyear , iidate
    real(rkx) :: w1 , w2
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'solar_irradiance'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( isolconst == 1 ) then
      solar_irradiance = 1367.0_rkx
    else
      if ( scenario(1:3) == 'SSP' ) then
        if ( .not. associated(heppatsi) ) then
          call read_solarforcing( )
        end if
        iidate = (rcmtimer%year-1850)*12 + rcmtimer%month
        if ( iidate < 1 ) then
          iidate = mod(iidate,132)+1
        else if ( iidate > 5400) then
          iidate = 5400 - 132 + mod(iidate,132)
        end if
        solar_irradiance = tsifac*heppatsi(iidate)
      else
        if ( calday > dayspy/2.0_rkx ) then
          w2 = calday/dayspy-0.5_rkx
          w1 = 1.0_rkx-w2
          iyear = rcmtimer%year
        else
          w1 = 0.5_rkx-calday/dayspy
          w2 = 1.0_rkx-w1
          iyear = rcmtimer%year-1
        end if
        iidate = rcmtimer%year*10000+rcmtimer%month*100+rcmtimer%day
        if ( iidate > 20080630 ) then
          iyear = mod(rcmtimer%year,12)+1996
        end if
        if ( iidate < 16100101 ) then
          iyear = 1610 + mod(rcmtimer%year,12)
        end if
        solar_irradiance = tsifac*(w1*tsi(3,iyear)+w2*tsi(3,iyear+1))
      end if
    end if
    if ( itweak == 1 ) then
      if ( itweak_solar_irradiance == 1 ) then
        solar_irradiance = solar_irradiance + solar_tweak
      end if
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end function solar_irradiance

  subroutine read_solarforcing( )
    implicit none
    character(len=256) :: heppafile
    integer(ik4) :: iret , ncid , idimid , ntime , ivar

    if ( myid == italk ) then
      write(stdout,*) 'Read solar forcing total irradiance data...'
    end if
    heppafile = trim(inpglob)//pthsep//'CMIP6'//pthsep//'SOLAR'//pthsep// &
      'solarforcing-ref-mon_input4MIPs_solar_CMIP_SOLARIS-HEPPA'// &
      '-3-2_gn_185001-229912.nc'
    iret = nf90_open(heppafile,nf90_nowrite,ncid)
    if ( iret /= nf90_noerr ) then
      write (stderr, *) trim(heppafile) , ': ', nf90_strerror(iret)
      call fatal(__FILE__,__LINE__,'CANNOT OPEN SOLARFORCING FILE')
    end if
    iret = nf90_inq_dimid(ncid,'time',idimid)
    if ( iret /= nf90_noerr ) then
      write (stderr, *) trim(heppafile) , ': ', nf90_strerror(iret)
      call fatal(__FILE__,__LINE__,'CANNOT FIND TIME DIMENSION')
    end if
    iret = nf90_inquire_dimension(ncid,idimid,len=ntime)
    if ( iret /= nf90_noerr ) then
      write (stderr, *) trim(heppafile) , ': ', nf90_strerror(iret)
      call fatal(__FILE__,__LINE__,'CANNOT READ TIME DIMENSION LENGHT')
    end if
    iret = nf90_inq_varid(ncid,'tsi',ivar)
    if ( iret /= nf90_noerr ) then
      write (stderr, *) trim(heppafile) , ': ', nf90_strerror(iret)
      call fatal(__FILE__,__LINE__,'VARIABLE TSI NOT FOUND')
    end if
    call getmem1d(heppatsi,1,ntime,'sun: heppatsi')
    iret = nf90_get_var(ncid,ivar,heppatsi)
    if ( iret /= nf90_noerr ) then
      write (stderr, *) trim(heppafile) , ': ', nf90_strerror(iret)
      call fatal(__FILE__,__LINE__,'CANNOT READ VARIABLE TSI')
    end if
    iret = nf90_close(ncid)
    if ( iret /= nf90_noerr ) then
      write (stderr, *) trim(heppafile) , ': ', nf90_strerror(iret)
      call fatal(__FILE__,__LINE__,'CANNOT CLOSE SOLARFORCING FILE')
    end if
  end subroutine read_solarforcing

end module mod_sun
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
