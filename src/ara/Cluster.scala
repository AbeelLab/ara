package ara

object Cluster {
  class Cluster(s: String) {

    def hasReference: Boolean = s match {
      case "L1-L2-L3-L4" => true
      case "L2-L3-L4" => true
      case "L4" => true
      case "L4.2" => true
      case "L4.2.2" => true
      case "L4.2.2.1" => true
      case "L4.2.2.1.1" => true
      case "L4.2.2.1.1.2" => true
      case "L4.2.2.1.1.2.2" => true
      case "L4.2.2.1.1.2.2.1" => true
      case "L4.2.2.1.1.2.2.1.1" => true
      case "L4.2.2.1.1.2.2.1.1.1" => true
      case "L4.2.2.1.1.2.2.1.1.1.1" => true
      case "L4.2.2.1.1.2.2.1.1.1.1.1" => true
      case _ => false
    }
    
  }
  
  implicit def seqtoBool(s: String) = new Cluster(s)
}