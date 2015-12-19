package ara

/**
 * Lines in VCF file of samples mapped against minimized reference genome
 */

class DRVcfLine(line: String) {
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

  override def toString(): String = region + "\t" + chrPos + "\t" + ref + "/" + alt + "\t" + filter
}

object DRVcfLine {

  /** Only filter SNPs */
  def unapply(line: String): Option[String] = {
    val arr = line.mkString.split("\t")
    val ref = arr(3)
    val alt = arr(4)
    val nucleotides = Array[String]("A", "C", "T", "G")
    if (nucleotides.contains(ref) && nucleotides.contains(alt)) {
      val info = arr(7).split(";")
      val bc = info(5)
      val qp = info(6)
      val ac = info(11)
      Some(line)
    }
    else None
  }

}