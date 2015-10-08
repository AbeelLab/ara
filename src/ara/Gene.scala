package ara

/**
 * Unapply to read Genes from gff-files 
 */

class Gene(val tag: String, val feature: String, val start: Int, val end: Int, val dir: String) {
  override def toString(): String = "Gene(" + tag + ", " + start + ", " + end + ", " + dir + ")"
}

object Gene {
  def unapply(s: String): Option[(String, String, Int, Int, String)] = {
    val arr = s.split("\t")
    val feature = arr(2)
    val start = arr(3).toInt
    val end = arr(4).toInt
    val direction = arr(6)
    val locus_tag = arr(8).split(";").filter(_.startsWith("locus_tag")).mkString.split("=")(1)
    Some(locus_tag.substring(1, locus_tag.size - 1), feature, start, end, direction)
  }
}
