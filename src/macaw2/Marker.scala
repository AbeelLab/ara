package macaw2

abstract class Marker

object Marker {
  def unapply(s: String): Option[(String, Int, Boolean)] = {
    val line = s.split("\t")
    val mutInfo = line(0)
    val count = line(1).toInt
    val presence = if (line(3) == "present") true else false
    Some(mutInfo, count, presence)
  }
}

case class DrugMarker(val mutInfo: String, val count: Int, val isPresent: Boolean) extends Marker {
  override def toString(): String = mutInfo + "\t" + count + "\t" + (if (isPresent) "present" else "absent")
  val arr = mutInfo.split("_")
  val snp = arr(0)
  val coordinate = snp.drop(1).dropRight(1).toInt
  val alt = snp.drop(snp.length - 1)
  val locus = arr(1)
  val sOrNs = arr(2)
  val markerType = arr(3)
  val drugs = arr(4)
}

case class ClusterMarker(val mutInfo: String, val count: Int, val isPresent: Boolean) extends Marker {
  val lineage = mutInfo.split("_")(2)
  override def toString(): String = mutInfo + "\t" + count + "\t" + (if (isPresent) "present" else "absent") + "\t" + lineage
}

