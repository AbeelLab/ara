package ara.dr

/**
 * Lines in VCF file of samples mapped against minimized reference genome
 */

class DRVcfLine(val line: String) {
  val arr = line.split("\t")
  val regionID = arr(0)
  val pos = arr(1).toInt
  val ref = arr(3)
  val alt = arr(4)
  val filter = arr(6)
  val regionArr = regionID.split("_")
  val region = regionArr(4)
  val regionStart = regionArr.last.split("-")(0).toInt
  val chrPos = regionStart + pos - 1
  val loci = region.split("/")

  lazy val info = arr(7).split(";")
  lazy val bc = info(5)
  lazy val qp = info(6)
  lazy val ac = info(11)

  val nucleotides = Array[String]("A", "C", "T", "G")
  def isSNP(): Boolean = nucleotides.contains(ref) && nucleotides.contains(alt)  
  
  def ambiguous: Boolean = filter == "Amb"
  def pass: Boolean = filter == "PASS"
  
  override def toString(): String = region + "\t" + chrPos + "\t" + ref + "/" + alt + "\t" + filter
  
  
}

object DRVcfLine {

  def apply(line: String): DRVcfLine = {
    new DRVcfLine(line)
  }
  
  /** Only filter SNPs */
  /*def unapply(line: String): Option[String] = {
    val arr = line.mkString.split("\t")
    val ref = arr(3)
    val alt = arr(4)
    if (nucleotides.contains(ref) && nucleotides.contains(alt)) {
      Some(line)
    } else None
  }*/

}