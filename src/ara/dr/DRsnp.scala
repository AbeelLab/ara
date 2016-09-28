package ara.dr

/**
 * Read drug resistance SNP from Coll2015
 */
class DRsnp(val drug: String, val locus: String, val locusTag: String, val cp: Int, val r: String, val gp: Int, val a: String, val codonNum: String, val aaChange: String) extends Ordered[DRsnp] {
  
  override def toString(): String = ">" + r + gp + a + "_" + locus + "_" + codonNum + "_" + aaChange + "_resistance_" + drug
  def complementToString(drugs: String): String = ">" + r + cp + r + "_susceptibility_" + drugs
  def compare(that: DRsnp): Int = this.cp compare that.cp
}

object DRsnp {
  def unapply(s: String): Option[(String, String, String, Int, String, Int, String, String, String)] = {
    val sArr = s.mkString.split("\t")
    val chrPos = sArr(2)
    if (!chrPos.contains("/") && !chrPos.equals("-")) {
      val drug = sArr(0)
      val locus = if (sArr(1).contains("_")) sArr(1).split("_").mkString("-") else sArr(1)
      val locusTag = sArr(8)
      val genePos = sArr(3)
      val nChange = sArr(4)
      val ncArr = nChange.split("/")
      val r = ncArr(0)
      val a = ncArr(1)
      val codonNumber = sArr(5)
      val aaChange = sArr(7)
      val nucleotides = Array[String]("A", "C", "T", "G")
      if (nucleotides.contains(r) && nucleotides.contains(a))
        Some((drug, locus, locusTag, chrPos.toInt, r, genePos.toInt, a, codonNumber, aaChange))
      else None
    } else None
  }

  class Line(s: String) {
    def isSNP: Boolean = s match {
      case DRsnp(d, l, lt, cp, r, gp, a, cn, aac) => true
      case _ => false
    }
  }

  implicit def seqtoBool(s: String) = new Line(s)
}

